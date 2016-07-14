#!~/Downloads/julia-0.4.5/usr/bin/julia
using JuMP
using CPLEX


function DANminus1proxy(N, Dn, D, Gn, Gtype, G, Plim, RD, RU, DT, UT, su_cost, cost, beta, X, fmax, STfmax, L, GenSchOut, GenForceOut, BranchSchOut, BranchForceOut, MarketClearing, Ploadfc, HHourCap, days)

    #DEFINE time constants
    hpd = 24
    dpw = 7
    wpy = 52
    hpy = hpd * dpw * wpy
    dpy = dpw * wpy
    Dt = 60                 #minutes
    dt = 15                 #minutes

    #CREATE N-1 CONTINGENCY LIST
    C = L
    a = ones(Float64, L, C) - eye(L, C)

    #create matrix of generator availability IN DAYS
    aGen = ones(Int, G, dpy)

    for g=1:G
        for week in GenSchOut[g]    #in weeks
            sch_start = dpw * (week - 1) + 1
            sch_end = min(dpy, sch_start + dpw)
            aGen[g, sch_start:sch_end] = 0
        end
        for (start, duration) in GenForceOut[g] #in hours
            f_start = Int(ceil(start/24)) + 1   #make it out of operation for the next day (day of failure plus one)
            f_end = min(hpy, Int(ceil((start+duration)/24)))
            if f_start <= f_end     #check if it's not fixed before the next day
                aGen[g, f_start:f_end] = 0
            end
        end
    end


    #next create matrix of branch availability IN DAYS
    aBranch = ones(Int, L, dpy)

    for l=1:L
        for (sch_start, sch_duration) in BranchSchOut[l]    #in days
            aBranch[l, sch_start:sch_start + sch_duration - 1] = 0
        end
        for (fcd_start, fcd_duration) in BranchForceOut[l]  #in hours
            d_start = Int(ceil(fcd_start/24)) + 1   #make it out of operation for the next day (day of failure plus one)
            d_end = min(hpy, Int(ceil((fcd_start+fcd_duration)/24)))
            if d_start <= d_end     #check if it's not fixed before the next day
                aBranch[l, d_start:d_end] = 0
            end
        end
    end


    #next create matrices of market clearing outcome
    uMC = zeros(Int, G, hpy)
    vMC = zeros(Int, G, hpy)
    PMC = zeros(Float64, G, hpy)
    for g=1:G
        for (h, v, Pgen) in MarketClearing[g]
            uMC[g,h] = 1
            vMC[g,h] = v
            PMC[g,h] = Pgen
        end
    end


    #DA DECISION and COST
    #create empty array indexed by generator, to be filled with tuples (hour, power output)
    DAdecision = Array(Any, G)
    for g=1:G
        DAdecision[g] = Tuple{Int, Float64}[]
    end
    #create array of costs
    DAcost = zeros(Float64, dpy)


    #initialize gen status and time of 'must keep status'
    aDA = zeros(Int, G)
    KTDA = zeros(Int, G)

    #set startup and shutdown ramp rates equal to Plim
    RSD = Dict{AbstractString,Float64}()
    RSU = Dict{AbstractString,Float64}()
    for k in keys(RD)
        RSD[k] = Plim[(k,5)]
        RSU[k] = RSD[k]
    end


    for day in days

        m = Model(solver = CplexSolver(CPX_PARAM_MIPDISPLAY=1, CPX_PARAM_EPGAP=0.5))


        @variable(m, uDA[1:G,1:hpd], Bin)
        @variable(m, 0 <= vDA[1:G,1:hpd] <= 1)
        @variable(m, 0 <= wDA[1:G,1:hpd] <= 1)
        @variable(m, vSO[1:G,1:hpd] >= 0)
        @variable(m, Pup0[1:G,1:hpd] >= 0)
        @variable(m, Pdn0[1:G,1:hpd] >= 0)
        @variable(m, Pupc[1:G,1:C,1:hpd] >= 0)
        @variable(m, Pdnc[1:G,1:C,1:hpd] >= 0)
        @variable(m, f0[1:L,1:hpd])
        @variable(m, th0[1:N,1:hpd])
        @variable(m, STfc[1:L,1:C,1:hpd])
        @variable(m, STthc[1:N,1:C,1:hpd])
        @variable(m, fc[1:L,1:C,1:hpd])
        @variable(m, thc[1:N,1:C,1:hpd])


        @objective(m, Min, sum{sum{su_cost[Gtype[g]] * vSO[g,h] + cost[(Gtype[g],4)] * (Pup0[g,h]) * 100, g=1:G}, h = 1:hpd} )

        #Generator availability and cycling
        @constraint(m, keep[g=1:G, h = 1:KTDA[g]], uDA[g,h] == aDA[g])
        @constraint(m, updown1[g=1:G], vDA[g,1] - wDA[g,1] == uDA[g,1] - aDA[g])
        @constraint(m, updown[g=1:G, h = 2:hpd], vDA[g,h] - wDA[g,h] == uDA[g,h] - uDA[g,h-1] )
        @constraint(m, MinUPTime[g=1:G, h = UT[Gtype[g]]:hpd], sum{vDA[g,s], s=h-UT[Gtype[g]]+1:h} <= uDA[g,h] )
        @constraint(m, MinDNTime[g=1:G, h = DT[Gtype[g]]:hpd], sum{wDA[g,s], s=h-DT[Gtype[g]]+1:h} <= 1 - uDA[g,h] )
        @constraint(m, TSOpay[g=1:G, h = 1:hpd], vSO[g,h] >= vDA[g,h] - vMC[g,h] )


        #Pre-contingency state
        @constraint(m, Gen0Min[g=1:G,h=1:hpd], Pup0[g,h] >= aGen[g] * uDA[g,h] * (Plim[(Gtype[g],1)] - PMC[g,hpd*(day-1)+h]))
        @constraint(m, noHydroGen0Max[g=1:G,h=1:hpd; Gtype[g] != "U50"], Pup0[g,h] <= aGen[g] * uDA[g,h] * (Plim[(Gtype[g],5)] - PMC[g,hpd*(day-1)+h]))
        @constraint(m, HydroGen0Max[g=1:G,h=1:hpd; Gtype[g] == "U50"], Pup0[g,h] <= aGen[g] * uDA[g,h] * (HHourCap[hpd*(day-1)+h] - PMC[g,hpd*(day-1)+h]))
        @constraint(m, Gen0Minbis[g=1:G,h=1:hpd], Pdn0[g,h] <= aGen[g] * uDA[g,h] * uMC[g,hpd*(day-1)+h] * (PMC[g,hpd*(day-1)+h] - Plim[(Gtype[g],1)]))


        @constraint(m, Power0[n=1:N, h=1:hpd], 
            sum{uDA[g,h] * PMC[g,hpd*(day-1)+h] + Pup0[g,h] - Pdn0[g,h], g in Gn[n]} - 
            sum{beta[n,l] * f0[l,h], l=1:L} == 
            sum{Ploadfc[d,hpd*(day-1)+h], d in Dn[n]}
            )
        @constraint(m, RampDownInterPeriod[g=1:G, h=2:hpd], 
            (uDA[g,h-1] * PMC[g,hpd*(day-1)+h-1] + Pup0[g,h-1] - Pdn0[g,h-1]) - 
            (uDA[g,h] * PMC[g,hpd*(day-1)+h] + Pup0[g,h] - Pdn0[g,h]) <= 
            Dt * RD[Gtype[g]] * uDA[g,h] + RSD[Gtype[g]] * wDA[g,h]
            )
        @constraint(m, RampUpInterPeriod[g=1:G, h=2:hpd], 
            (uDA[g,h] * PMC[g,hpd*(day-1)+h] + Pup0[g,h] - Pdn0[g,h]) - 
            (uDA[g,h-1] * PMC[g,hpd*(day-1)+h-1] + Pup0[g,h-1] - Pdn0[g,h-1]) <= 
            Dt * RU[Gtype[g]] * uDA[g,h-1] + RSU[Gtype[g]] * vDA[g,h]
            )

        @constraint(m, Flow0[l=1:L, h=1:hpd], f0[l,h] - aBranch[l,day] * (1/X[l]) * sum{beta[n,l] * th0[n,h], n=1:N} == 0 )
        @constraint(m, Flow0PosLim[l=1:L, h=1:hpd], f0[l,h] <= fmax[l] )
        @constraint(m, Flow0NegLim[l=1:L, h=1:hpd], -f0[l,h] <= fmax[l] )


        @constraint(m, Reference0[h=1:hpd], th0[1,h] == 0 )
        @constraint(m, ReferenceSTc[c=1:C, h=1:hpd], STthc[1,c,h] == 0 )
        @constraint(m, Referencec[c=1:C, h=1:hpd], thc[1,c,h] == 0 )


        #Short-term post-contingency state
        @constraint(m, PowerSTc[n=1:N,c=1:C,h=1:hpd], 
            sum{uDA[g,h] * PMC[g,hpd*(day-1)+h] + Pup0[g,h] - Pdn0[g,h], g in Gn[n]} - 
            sum{beta[n,l] * STfc[l,c,h], l=1:L} == 
            sum{Ploadfc[d,hpd*(day-1)+h], d in Dn[n]}
            )
        @constraint(m, FlowSTc[l=1:L,c=1:C,h=1:hpd], STfc[l,c,h] - aBranch[l,day] * a[l,c] * (1/X[l]) * sum{beta[n,l] * STthc[n,c,h], n=1:N} == 0 )

        @constraint(m, FlowSTcPosLim[l=1:L,c=1:C,h=1:hpd], STfc[l,c,h] <= STfmax[l] )
        @constraint(m, FlowSTcNegLim[l=1:L,c=1:C,h=1:hpd], -STfc[l,c,h] <= STfmax[l] )


        #Post-contingency state (after the successful application of corrective control)
        @constraint(m, Powerc[n=1:N,c=1:C,h=1:hpd], 
            sum{uDA[g,h] * PMC[g,hpd*(day-1)+h] + Pup0[g,h] - Pdn0[g,h] + Pupc[g,c,h] - Pdnc[g,c,h], g in Gn[n]} - 
            sum{beta[n,l] * fc[l,c,h], l=1:L} == 
            sum{Ploadfc[d,hpd*(day-1)+h], d in Dn[n]}
            )
        @constraint(m, Flowc[l=1:L,c=1:C,h=1:hpd], fc[l,c,h] - aBranch[l,day] * a[l,c] * (1/X[l]) * sum{beta[n,l] * thc[n,c,h], n=1:N} == 0 )

        @constraint(m, FlowcPosLim[l=1:L,c=1:C,h=1:hpd], fc[l,c,h] <= fmax[l] )
        @constraint(m, FlowcNegLim[l=1:L,c=1:C,h=1:hpd], -fc[l,c,h] <= fmax[l] )
     
        @constraint(m, GencMin[g=1:G,c=1:C,h=1:hpd], Pup0[g,h] - Pdn0[g,h] + Pupc[g,c,h] >= aGen[g] * uDA[g,h] * (Plim[(Gtype[g],1)] - PMC[g,hpd*(day-1)+h]))
        @constraint(m, noHydroGencMax[g=1:G,c=1:C,h=1:hpd; Gtype[g] != "U50"], Pup0[g,h] - Pdn0[g,h] + Pupc[g,c,h] <= aGen[g] * uDA[g,h] * (Plim[(Gtype[g],5)] - PMC[g,hpd*(day-1)+h]))
        @constraint(m, HydroGencMax[g=1:G,c=1:C,h=1:hpd; Gtype[g] == "U50"], Pup0[g,h] - Pdn0[g,h] + Pupc[g,c,h] <= aGen[g] * uDA[g,h] * (HHourCap[hpd*(day-1)+h] - PMC[g,hpd*(day-1)+h]))
        @constraint(m, GencMinbis[g=1:G,c=1:C,h=1:hpd], Pdn0[g,h] - Pup0[g,h] + Pdnc[g,c,h] <= aGen[g] * uDA[g,h] * uMC[g,hpd*(day-1)+h] * (PMC[g,hpd*(day-1)+h] - Plim[(Gtype[g],1)]))

        @constraint(m, RampDownc[g=1:G,c=1:C,h=1:hpd], Pdnc[g,c,h] <= dt * RD[Gtype[g]] )
        @constraint(m, RampUpc[g=1:G,c=1:C,h=1:hpd], Pupc[g,c,h] <= dt * RU[Gtype[g]] )

        status = solve(m)

        uDAstar = round(Int, getvalue(uDA))
        PDAstar = getvalue(Pup0 - Pdn0)

        for g=1:G
            for h=1:hpd
                if uDAstar[g,h] == 1
                    hour = hpd * (day - 1) + h
                    push!(DAdecision[g], (hour, PMC[g,h] + PDAstar[g,h]))
                end
            end
        end   

        DAcost[day] = getobjectivevalue(m) 

        println("Solved in ", getsolvetime(m::Model))

    end

    return DAdecision, DAcost

end

