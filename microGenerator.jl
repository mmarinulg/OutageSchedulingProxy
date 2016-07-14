#!~/Downloads/julia-0.4.5/usr/bin/julia
tic()
using JuMP
using Cbc
using Gurobi
using CPLEX
using Distributions

T = []
push!(T, toq())
println("Task ", ARGS[2], ": Loading modules took ", T[end])


#DEFINE time constants
hpd = 24
dpw = 7
wpy = 52
hpy = hpd * dpw * wpy
dpy = dpw * wpy
hpw = hpd * dpw
Dt = 60                 #minutes


tic()
include("./readData.jl")
N = readBusData(string(ARGS[1]))
Dn, D = readLoadData(string(ARGS[1]), N)
Gn, Gtype, G = readGeneratorData(string(ARGS[1]), N)
Plim, RD, RU, DT, UT, su_cost, cost, GenOutRate, GenOutDuration, GenMaintPeriods = readGenLookupTable(string(ARGS[1]))
beta, X, fmax, STfmax, L, BranchOutRate, BranchOutDuration = readBranchData(string(ARGS[1]), N)
PCTWeeklyPeakLoad, PCTDailyPeakLoad, PCTHourlyPeakLoad, PCTBusLoad = readPeakLoadData(string(ARGS[1]))
push!(T, toq())
println("Task ", ARGS[2], ": Reading data took ", T[end])

#First simulate the annual peak

srand(parse(Int, ARGS[2]))	#seed from argument line
U = Uniform(27, 33) 		#uniformily distributed between 27 and 33 pu
#AnnualPeakLoad = rand(U)

AnnualPeakLoad = 30

#########################
#GENERATOR SCHEDULED OUTAGES (IN WEEKS)

tic()
#first create outage array of arrays 
GenSchOut = Array(Any, G)

#A) GENERATOR SCHEDULED OUTAGES (eHighWay methodology)

#create array of tuples (generator, capacity) and sort it by decreasing capacity
Gcap = [tuple(g, Plim[(Gtype[g],4)]) for g=1:G]
orderedGcap = sort(collect(Gcap), by=x->x[2], rev=true)

#initialize virtual residual load (using a weekly basis)
VRLoad = (.01) * PCTWeeklyPeakLoad * AnnualPeakLoad

#main loop
for (g,cap) in orderedGcap
	#println(g, ":", nb_maint[Gtype[g]])
	#Compute number of maintenance periods (of 1 week length)
	n = GenMaintPeriods[Gtype[g]]
	#create array that will be filled with the outages (in weeks)
	GenSchOutTemp = Int[]
	#create a temporary copy of the virtual residual load and compute its maximum (cheat)
	tempVRLoad = VRLoad
	max = maximum(tempVRLoad)
	#loop over maintenance periods
	for m = 1:n
		#pick the week with the lowest virtual residual load
		w = indmin(tempVRLoad)
		#schedule the maintenance
		push!(GenSchOutTemp, w)
		#make sure the same week doesn't get picked again
		tempVRLoad[w] = tempVRLoad[w] + max
	end
	GenSchOut[g] = sort(GenSchOutTemp)
	#now increase VRLoad by the capacity of g in the scheduled periods
	VRLoad[GenSchOut[g]] = VRLoad[GenSchOut[g]] + cap

end
push!(T, toq())
println("Task ", ARGS[2], ": Generator scheduled outages took ", T[end])

#GENERATOR FORCED OUTAGES (IN HOURS)

tic()
#first create outage array of arrays 
GenForceOut = Array(Any, G)

for g=1:G
	#Compute hourly outage rate (outages per hour)
	hGenOutRate = GenOutRate[Gtype[g]] / (hpd * 365)
	#Compute hourly outage probability
	hGenOutProb = hGenOutRate * exp(-hGenOutRate)
	#Compute outage expected duration (in hours)
	ExpGenOutDur = GenOutDuration[Gtype[g]]
	#create array that will be filled with the tuples 'outage time, duration' in hours
	GenForceOut[g] = Tuple{Int,Int}[]
	#start sweeping all hours in the year
	hour = 1
	while hour <= hpy
		day = Int(ceil(hour/24))
		week = Int(ceil(day/7))
		#check if generator already in scheduled maintenance for current week
		if week in GenSchOut[g]
			hour = hour + hpw	#go to next week
			continue
		end
		#simulate forced outage
		if rand(Uniform()) <= hGenOutProb
			#simulate outage duration
			GenOutDur = round(Int, rand(Exponential(ExpGenOutDur)))
			#add tuple
			push!(GenForceOut[g], (hour, GenOutDur))
			hour = hour + GenOutDur
			continue
		end
		hour = hour + 1	
	end
end
push!(T, toq())
println("Task ", ARGS[2], ": Generator forced outages took ", T[end])



#BRANCH FORCED OUTAGES (IN HOURS)

tic()
#first create outage array of arrays 
BranchForceOut = Array(Any, L)
#Compute hourly outage rate (outages per hour)
hBranchOutRate = BranchOutRate / hpy

for l=1:L
	#Compute hourly outage probability
	hBranchOutProb = hBranchOutRate[l] * exp(-hBranchOutRate[l])
	#create array that will be filled with the tuples 'outage time, duration' in hours
	BranchForceOut[l] = Tuple{Int, Int}[]
	#start sweeping all hours in the year
	hour = 1
	while hour <= hpy
		day = Int(ceil(hour/24))
		#assume scheduled maintenance is not decided yet		
		#simulate forced outage
		if rand(Uniform()) <= hBranchOutProb
			#simulate outage duration
			BranchOutDur = round(Int, rand(Normal(BranchOutDuration[l], 1)))
			#add tuple
			push!(BranchForceOut[l], (hour, BranchOutDur))
			hour = hour + BranchOutDur
			continue
		end
		hour = hour + 1	

	end
end
push!(T, toq())
println("Task ", ARGS[2], ": Branch forced outages took ", T[end])


####################
#DEMAND FORECAST
tic()

Ploadfc = []
PCTHourlyLoad = [] 	#for hydro (see below)

for hour = 1:hpy
	#println("Computing demand realization for hour ", hour)
	day = Int(ceil(hour/24))
	week = Int(ceil(day/7))
	hour_of_day = mod(hour-1,hpd) + 1
	day_of_week = mod(day-1,dpw) + 1
	day_type = day_of_week <= 5 ? 1:2		#day_type: 1 = weekday, 2 = weekend
	if week <= 8 || week >= 44
  		week_type = 0
	elseif week >= 18 && week <= 30
  		week_type = 2
	else
  		week_type = 4
	end 									#week_type: 0 = winter, 2 = summer, 4 = spring or fall
	push!(PCTHourlyLoad, (.01^3) * PCTHourlyPeakLoad[hour_of_day, week_type + day_type] * PCTDailyPeakLoad[day_of_week] * PCTWeeklyPeakLoad[week])
	temp = (.01) * PCTHourlyLoad[end] * AnnualPeakLoad * PCTBusLoad
	append!(Ploadfc, temp)
end
Ploadfc = reshape(Ploadfc, D, hpy)
push!(T, toq())

println("Task ", ARGS[2], ": Demand forecast took ", T[end])


#################
#DEMAND REALISATION

tic()

Pload = []

#Laurine's method
sigma_glob = .03
sigma_loc = 0.1
BusLoad = .01 * PCTBusLoad
sum_squares = sum(BusLoad .* BusLoad)

alpha_d = Float16(sqrt((sigma_glob^2 - sum_squares * sigma_loc^2) / (1 - sum_squares)))
beta_d = Float16(sqrt((sigma_glob^2 - sigma_loc^2) / (sum_squares - 1)))

N_alpha = Normal(0, alpha_d)
N_beta = Normal(0, beta_d)

for hour = 1:hpy
	eps_alpha = rand(N_alpha)
	eps_beta = rand(N_beta, D)
	temp = (1 + eps_alpha + eps_beta) .* Ploadfc[:,hour]
	append!(Pload, temp)
end
Pload = reshape(Pload, D, hpy)
push!(T, toq())
println("Task ", ARGS[2], ": Demand realization took ", T[end])


#################
#CAPACITY OF HYDRO

tic()
#first compute capacity and energy for different quarters. This is taken from the footnote in Table 9 of RTS-96 specs. Kinda hard-wired for now...
HQuarterCap = [1.0 1.0 0.9 0.9] * Plim[("U50",2)]	#hydro max is in segment k = 2
HQuarterEn = [.35 .35 .1 .2] * 200000 / 100 		#in pu

HHourCap = zeros(Float64, hpy)

#initialize first and last hour of the quarter
hpq = round(Int, hpy/4)
first = 1
last = hpq
for quarter = 1:4
	#create array of tuples (hour, load) for the quarter
	HourLoad = [tuple(h, PCTHourlyLoad[h]) for h=first:last]
	#sort the array by load in decreasing order
	ordHourLoad = sort(collect(HourLoad), by=x->x[2], rev=true)
	#compute capacity and available energy for the quarter
	quarterCap = HQuarterCap[quarter]
	quarterEn = HQuarterEn[quarter]
	for (h, load) in ordHourLoad
		#compute capacity to be assigned as the minimum between the capacity of the quarter and the energy available
		cap = min(quarterCap, quarterEn)
		#assign capacity
		HHourCap[h] = cap
		#reduce available energy
		quarterEn = quarterEn - cap
		if quarterEn <= 0 break end
	end
	#compute first and last hour of the next quarter
	first = first + hpq
	last = last + hpq

end

push!(T, toq())
println("Task ", ARGS[2], ": Hydro capacity took ", T[end])


#########################
#MARKET CLEARING OUTCOME

#first create empty array indexed by generator, to be filled with tuples (hour, startup variable, power output)
MarketClearing = Array(Any, G)
for g=1:G
	MarketClearing[g] = Tuple{Int, Int, Float64}[]
end


#initialize gen status and time of 'must keep status'
aMC = zeros(Int, G)
KTMC = zeros(Int, G)

#set startup and shutdown ramp rates equal to Plim
RSD = Dict{AbstractString,Float64}()
RSU = Dict{AbstractString,Float64}()
for k in keys(RD)
	RSD[k] = Plim[(k,5)]
	RSU[k] = RSD[k]
end

#create vector of generator availability
aGen = ones(Int, G)

#for day=1:dpy
for day=1:2

    tic()

	#fill vector of generator availability
	
	week = Int(ceil(day/7))
	hour = hpd * (day - 1)	#compute last hour in the previous day, to check for forced outages
    for g=1:G
    	if week in GenSchOut[g]
    		aGen[g] = 0
    		continue
    	end
    	for (start, duration) in GenForceOut[g]
    		if start > hour break		#if the current outage starts after the considered hour, no need to look at more outages (since GenForceOut is ordered)
    		elseif hour <= start + duration 	#if the generator is out at the considered hour, it is out for the entire day (conservative approach)
    			aGen[g] = 0
    		end
    	end
    end
    
	
	#m = Model(solver = GurobiSolver(OutputFlag=0))
    #m = Model(solver = CbcSolver(logLevel=0, threads=16))
	m = Model(solver = CplexSolver(CPX_PARAM_MIPDISPLAY=0, CPX_PARAM_EPGAP=0.001))
	#m = Model(solver = CplexSolver())

	@variable(m, PMC[1:G,1:hpd])
	@variable(m, PMCk[1:G,1:hpd,1:4] >= 0)
	@variable(m, uMC[1:G,1:hpd], Bin)
	@variable(m, 0 <= vMC[1:G,1:hpd] <= 1)
	@variable(m, 0 <= wMC[1:G,1:hpd] <= 1)

	@objective(m, Min, sum{sum{su_cost[Gtype[g]] * vMC[g,h] + sum{cost[(Gtype[g],k)] * PMCk[g,h,k] * 100, k=1:4}, g=1:G}, h = 1:hpd} )

    #Generator availability and cycling
	@constraint(m, keep[g=1:G, h = 1:KTMC[g]], uMC[g,h] == aMC[g])
	@constraint(m, updown1[g=1:G], vMC[g,1] - wMC[g,1] == uMC[g,1] - aMC[g])
	@constraint(m, updown[g=1:G, h = 2:hpd], vMC[g,h] - wMC[g,h] == uMC[g,h] - uMC[g,h-1] )
	@constraint(m, MinUPTime[g=1:G, h = UT[Gtype[g]]:hpd], sum{vMC[g,s], s=h-UT[Gtype[g]]+1:h} <= uMC[g,h] )
	@constraint(m, MinDNTime[g=1:G, h = DT[Gtype[g]]:hpd], sum{wMC[g,s], s=h-DT[Gtype[g]]+1:h} <= 1 - uMC[g,h] )

	#Network balance and limits
	@constraint(m, NetworkBalance[h=1:hpd], sum{aGen[g] * PMC[g,h], g=1:G} ==
		sum{Ploadfc[d, hpd * (day - 1) + h], d=1:D} )

	@constraint(m, AllMaxPLimk1[g=1:G, h=1:hpd], PMCk[g,h,1] == uMC[g,h] * Plim[(Gtype[g],1)] )
	@constraint(m, HydroMaxPLimk2[g=1:G, h=1:hpd; Gtype[g] == "U50"], PMCk[g,h,2] <= uMC[g,h] * HHourCap[hpd*(day-1)+h] ) 
	@constraint(m, noHydroMaxPLimk2[g=1:G, h=1:hpd; Gtype[g] != "U50"], PMCk[g,h,2] <= uMC[g,h] * Plim[(Gtype[g],2)] ) 
	@constraint(m, AllMaxPLimk34[g=1:G, h=1:hpd, k=3:4], PMCk[g,h,k] <= uMC[g,h] * Plim[(Gtype[g],k)] ) 
	@constraint(m, sumPk[g=1:G, h = 1:hpd], PMC[g,h] == sum{PMCk[g,h,k], k=1:4})

	@constraint(m, RampDown[g=1:G, h=2:hpd], PMC[g,h-1] - PMC[g,h] <= Dt * RD[Gtype[g]] * uMC[g,h] + RSD[Gtype[g]] * wMC[g,h])
	@constraint(m, RampUp[g=1:G, h=2:hpd], PMC[g,h] - PMC[g,h-1] <= Dt * RU[Gtype[g]] * uMC[g,h-1] + RSU[Gtype[g]] * vMC[g,h])

	status = solve(m)

	uMCstar = round(Int, getvalue(uMC))
	vMCstar = round(Int, getvalue(vMC))
	PMCstar = getvalue(PMC)

	for g=1:G
		for h=1:hpd
			if uMCstar[g,h] == 1
				hour = hpd * (day - 1) + h
				push!(MarketClearing[g], (hour, vMCstar[g,h], PMCstar[g,h]))
			end
		end
	end

    push!(T, toq())

    println("Task ", ARGS[2], ": Computing market outcome for day ", day, " took ", T[end])
    println("Task ", ARGS[2], ": Objective value = ", getobjectivevalue(m))
    #println("Task ", ARGS[2], ": Power output = ", getvalue(PMC))

end


#Write data to file
tic()

outfile = string("./micro-scenarios/micro",ARGS[2],".",ARGS[1])
if(isfile(outfile)) rm(outfile) end

open(outfile,"a") do x
	write(x,"GENERATOR SCHEDULED OUTAGES\n")
	for g=1:G
		if isempty(GenSchOut[g]) continue end	#skip generators with no outages
		write(x, string(g), " ")	#write generator index	
		write(x, join(GenSchOut[g], " "))
		write(x, "\n")
	end
	write(x,"END\n\n")

	write(x,"GENERATOR FORCED OUTAGES\n")
	for g=1:G
		if isempty(GenForceOut[g]) continue	end 	#skip generators with no outages
		write(x, string(g), " ")	#write generator index	
		write(x, join(GenForceOut[g], " "))
		write(x, "\n")
	end	
	write(x,"END\n\n")

	write(x,"BRANCH FORCED OUTAGES\n")
	for l=1:L
		if isempty(BranchForceOut[l]) continue end 	#skip branches with no outages
		write(x, string(l), " ")	#write branch index	
		write(x, join(BranchForceOut[l], " "))
		write(x, "\n")
	end	
	write(x,"END\n\n")

	write(x,"DEMAND FORECAST\n")
	for d=1:D
		write(x, join(Ploadfc[d,:], " " ))
		write(x, "\n")
	end
	write(x,"END\n\n")

	write(x,"DEMAND REALIZATION\n")
	for d=1:D
		write(x, join(Pload[d,:], " " ))
		write(x, "\n")
	end
	write(x,"END\n\n")

	write(x,"HYDRO CAPACITY\n")
	write(x, join(HHourCap, " "))
	write(x, "\n")
	write(x,"END\n\n")

	write(x,"MARKET CLEARING OUTCOME\n")
	for g=1:G
		if isempty(MarketClearing[g]) continue	end 	#skip generators which will remain off
		write(x, string(g), " ")	#write generator index		
		write(x, join(MarketClearing[g], " "))
		write(x, "\n")
	end
	write(x,"END\n\n")
	
	
end

push!(T, toq())
println("Task ", ARGS[2], ": Writing output file took ", T[end])

println("Task ", ARGS[2], ": Everything took ", sum(T))

