#!~/Downloads/julia-0.4.5/usr/bin/julia
using JuMP
using CPLEX


function RTNminus1proxy(N, Dn, D, Gn, Gtype, G, Plim, RD, RU, cost, beta, X, fmax, STfmax, L, BranchOutRate, BranchSchOut, BranchForceOut, DAdecision, Pload, HHourCap, VoLL, hours)
    #DEFINE time constants
    hpd = 24
    dpw = 7
    wpy = 52
    hpy = hpd * dpw * wpy
    dpy = dpw * wpy
    Dt = 60                 #minutes
    dt = 15                 #minutes


    #Compute hourly outage rate (outages per hour)
    hBranchOutRate = BranchOutRate / hpy
    hAnyBranchOutProb = 0
    for l=1:L
        #Compute hourly outage probability
        hAnyBranchOutProb = hAnyBranchOutProb + hBranchOutRate[l] * exp(-hBranchOutRate[l])
    end
    println(hAnyBranchOutProb)

    #create matrix of branch availability IN HOURS
    uBranch = ones(Int, L, hpy)
    for l=1:L
        for (start, duration) in BranchSchOut[l]    #in days
            aBranch[l, start:start + duration] = 0
        end
        for (start, duration) in BranchForceOut[l]
            stop = min(hpy, start + duration)
            uBranch[l, start:stop] = 0
        end
    end

    #next create matrices of DA decisions
    uDA = zeros(Int, G, hpy)
    PDA = zeros(Float64, G, hpy)
    for g=1:G
        for (h, Pgen) in DAdecision[g]
            uDA[g,h] = 1
            PDA[g,h] = Pgen
        end
    end


    #CREATE N-1 CONTINGENCY LIST
    C = L
    a = ones(Float64, L, C) - eye(L, C)

    RTcost = zeros(Float64, hpy)

    #for h = 1:hpy
    for h in hours

        tic()
        #m = Model(solver = GurobiSolver(OutputFlag=0))
        #m = Model(solver = GLPKSolverLP())
        #m = Model(solver = ClpSolver())
        #m = Model()
        m = Model(solver = CplexSolver(CPX_PARAM_SIMDISPLAY=0))

        @variable(m, Pup0[1:G] >= 0)
        @variable(m, Pdn0[1:G] >= 0)
        @variable(m, Pupc[1:G,1:C] >= 0)
        @variable(m, Pdnc[1:G,1:C] >= 0)
        @variable(m, f0[1:L])
        @variable(m, th0[1:N])
        @variable(m, STfc[1:L,1:C])
        @variable(m, STthc[1:N,1:C])
        @variable(m, fc[1:L,1:C])
        @variable(m, thc[1:N,1:C])
        @variable(m, PloadSH[1:D] >= 0)

        for d=1:D
            setupperbound(PloadSH[d], 0.0)
        end

        @objective(m, Min, sum{cost[(Gtype[g],4)] * (Pup0[g] - Pdn0[g]) * 100, g=1:G} )

        #Pre-contingency state
        @constraint(m, Power0[n=1:N], 
            sum{PDA[g,h] + (Pup0[g] - Pdn0[g]), g in Gn[n]} - 
            sum{beta[n,l] * f0[l], l=1:L} == 
            sum{Pload[d,h] - PloadSH[d], d in Dn[n]}
            )
        @constraint(m, Flow0[l=1:L], f0[l] - uBranch[l,h] * (1/X[l]) * sum{beta[n,l] * th0[n], n=1:N} == 0 )

        @constraint(m, Flow0PosLim[l=1:L], f0[l] <= fmax[l] )
        @constraint(m, Flow0NegLim[l=1:L], -f0[l] <= fmax[l] )

        @constraint(m, GenMin0[g=1:G], PDA[g,h] - Pdn0[g] >= uDA[g,h] * Plim[(Gtype[g],1)])
        @constraint(m, GenMax0[g=1:G; Gtype[g] != "U50"], PDA[g,h] + Pup0[g] <= uDA[g,h] * Plim[(Gtype[g],5)])
        @constraint(m, GenMax0Hydro[g=1:G; Gtype[g] == "U50"], PDA[g,h] + Pup0[g] <= uDA[g,h] * HHourCap[h])

        @constraint(m, RampDown0[g=1:G], Pdn0[g] <= 0.5 * Dt * RD[Gtype[g]] )
        @constraint(m, RampUp0[g=1:G], Pup0[g] <= 0.5 * Dt * RU[Gtype[g]] )

        #Short-term post-contingency state
        @constraint(m, STPowerc[n=1:N,c=1:C], 
            sum{PDA[g,h] + (Pup0[g] - Pdn0[g]), g in Gn[n]} - 
            sum{beta[n,l] * STfc[l,c], l=1:L} == 
            sum{Pload[d,h] - PloadSH[d], d in Dn[n]}
            )
        @constraint(m, STFlowc[l=1:L,c=1:C], STfc[l,c] - uBranch[l,h] * a[l,c] * (1/X[l]) * sum{beta[n,l] * STthc[n,c], n=1:N} == 0 )

        @constraint(m, STFlowcPosLim[l=1:L,c=1:C], STfc[l,c] <= STfmax[l] )
        @constraint(m, STFlowcNegLim[l=1:L,c=1:C], -STfc[l,c] <= STfmax[l] )

        #Post-contingency state (after the successful application of corrective control)
        @constraint(m, Powerc[n=1:N,c=1:C], 
            sum{PDA[g,h] + (Pup0[g] - Pdn0[g]) + (Pupc[g,c] - Pdnc[g,c]), g in Gn[n]} - 
            sum{beta[n,l] * fc[l,c], l=1:L} == 
            sum{Pload[d,h] - PloadSH[d], d in Dn[n]}
            )
        @constraint(m, Flowc[l=1:L,c=1:C], fc[l,c] - uBranch[l,h] * a[l,c] * (1/X[l]) * sum{beta[n,l] * thc[n,c], n=1:N} == 0 )

        @constraint(m, FlowcPosLim[l=1:L,c=1:C], fc[l,c] <= fmax[l] )
        @constraint(m, FlowcNegLim[l=1:L,c=1:C], -fc[l,c] <= fmax[l] )

        @constraint(m, GenMinc[g=1:G,c=1:C], PDA[g,h] + Pup0[g] - Pdn0[g] - Pdnc[g,c] >= uDA[g,h] * Plim[(Gtype[g],1)])
        @constraint(m, GenMaxc[g=1:G,c=1:C; Gtype[g] != "U50"], PDA[g,h] + Pup0[g] - Pdn0[g] + Pupc[g,c] <= uDA[g,h] * Plim[(Gtype[g],5)])
        @constraint(m, GenMaxcHydro[g=1:G,c=1:C; Gtype[g] == "U50"], PDA[g,h] + Pup0[g] - Pdn0[g] + Pupc[g,c] <= uDA[g,h] * HHourCap[h])

        @constraint(m, RampDownc[g=1:G,c=1:C], Pdnc[g,c] <= dt * RD[Gtype[g]] )
        @constraint(m, RampUpc[g=1:G,c=1:C], Pupc[g,c] <= dt * RU[Gtype[g]] )


        push!(T, toq())
        println("Task ", ARGS[2], ": Setting up the SCOPF model for hour ", h, " took " , T[end])

        tic()
        status = solve(m; suppress_warnings=true)
        push!(T, toq())
        println("Task ", ARGS[2], ": Solving the SCOPF for hour ", h, " took " , T[end])

        if status == :Optimal
            RTcost[h] = getobjectivevalue(m)
            println("Task ", ARGS[2], ": RT cost = ", RTcost[h])
            continue
        end

        println("SCOPF not solved to optimality (infeasible?), activate preventive load shedding...")

        for d=1:D
            setupperbound(PloadSH[d], Pload[d,h])
        end
        @objective(m, Min, sum{cost[(Gtype[g],4)] * (Pup0[g] - Pdn0[g]) * 100, g=1:G} + sum{VoLL[d] * PloadSH[d] * 100, d=1:D} )
        status = solve(m; suppress_warnings=true)
        if status == :Optimal
            RTcost[h] = getobjectivevalue(m)
            println("Task ", ARGS[2], ": RT cost = ", RTcost[h])
            continue
        end

        println("SCOPF with load shedding not solved to optimality (infeasible?), falling back to OPF...")

        m = 0       #clear previous model

        m = Model(solver = CplexSolver(CPX_PARAM_SIMDISPLAY=0))

        @variable(m, Pup0[1:G] >= 0)
        @variable(m, Pdn0[1:G] >= 0)
        @variable(m, f0[1:L])
        @variable(m, th0[1:N])
        @variable(m, PloadSH[1:D] >= 0)

        @objective(m, Min, sum{cost[(Gtype[g],4)] * (Pup0[g] - Pdn0[g]) * 100, g=1:G} + sum{VoLL[d] * PloadSH[d] * 100, d=1:D} )

        #Pre-contingency state
        @constraint(m, Power0[n=1:N], 
            sum{PDA[g,h] + (Pup0[g] - Pdn0[g]), g in Gn[n]} - 
            sum{beta[n,l] * f0[l], l=1:L} == 
            sum{Pload[d,h] - PloadSH[d], d in Dn[n]}
            )
        @constraint(m, Flow0[l=1:L], f0[l] - uBranch[l,h] * (1/X[l]) * sum{beta[n,l] * th0[n], n=1:N} == 0 )

        @constraint(m, PloadMin[d=1:D], Pload[d,h] - PloadSH[d] >= 0 )

        @constraint(m, Flow0PosLim[l=1:L], f0[l] <= fmax[l] )
        @constraint(m, Flow0NegLim[l=1:L], -f0[l] <= fmax[l] )

        @constraint(m, GenMin0[g=1:G], PDA[g,h] - Pdn0[g] >= uDA[g,h] * Plim[(Gtype[g],1)])
        @constraint(m, GenMax0[g=1:G; Gtype[g] != "U50"], PDA[g,h] + Pup0[g] <= uDA[g,h] * Plim[(Gtype[g],5)])
        @constraint(m, GenMax0Hydro[g=1:G; Gtype[g] == "U50"], PDA[g,h] + Pup0[g] <= uDA[g,h] * HHourCap[h])

        @constraint(m, RampDown0[g=1:G], Pdn0[g] <= 0.5 * Dt * RD[Gtype[g]] )
        @constraint(m, RampUp0[g=1:G], Pup0[g] <= 0.5 * Dt * RU[Gtype[g]] )

        status = solve(m; suppress_warnings=true)

        if status == :Optimal
            PloadSHstar = getvalue(PloadSH)
            extra_cost = hAnyBranchOutProb * sum(VoLL .* (Pload[:,h] - PloadSHstar))
            RTcost[h] = getobjectivevalue(m) + extra_cost
            println("Task ", ARGS[2], ": RT cost = ", RTcost[h])
            continue
        end

        println("OPF with load shedding not solved to optimality (infeasible?), lose all system...")
        RTcost[h] =  sum(VoLL .* Pload[:,h])
        println("Task ", ARGS[2], ": RT cost = ", RTcost[h])


    end

    return RTcost

end

