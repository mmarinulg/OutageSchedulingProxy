#!~/Downloads/julia-0.4.5/usr/bin/julia

include("./readData.jl")
include("./DAN-1proxy.jl")
include("./RTN-1proxy.jl")

#DEFINE time constants
hpd = 24
dpw = 7
wpy = 52
hpy = hpd * dpw * wpy
dpy = dpw * wpy
Dt = 60                 #minutes
dt = 15                 #minutes


#READ ALL NECESSARY DATA

tic()
N = readBusData(string(ARGS[1]))
Dn, D = readLoadData(string(ARGS[1]), N)
Gn, Gtype, G = readGeneratorData(string(ARGS[1]), N)
Plim, RD, RU, DT, UT, su_cost, cost = readGenLookupTable(string(ARGS[1]))
beta, X, fmax, STfmax, L, BranchOutRate = readBranchData(string(ARGS[1]), N)
GenSchOut, GenForceOut = readGenAvailability(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]), G)
BranchForceOut = readBranchAvailability(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]), L)
MarketClearing = readMarketOutcome(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]), G)
Ploadfc = readDemandForecast(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]))
Pload = readDemandRealization(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]))
HHourCap = ReadHydroCapacity(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]))
VoLL = readVoLLData(string(ARGS[1]))
LineOutages, LineOutDuration = readMaintPolicyData("./maintenance_policy.csv")


#create outage array of arrays
BranchSchOut = Array{Any}(L)
for l = 1:L
	BranchSchOut[l] = Tuple{Int, Int}[]
end


for l in LineOutages
	#define outage spanning all year
	BranchSchOut[l] = [(1, dpy)]
	#compute associated DA decisions and costs
	DAdecision, DAcost = DANminus1proxy(N, Dn, D, Gn, Gtype, G, Plim, RD, RU, DT, UT, su_cost, cost, beta, X, fmax, STfmax, L, GenSchOut, GenForceOut, BranchSchOut, BranchForceOut, MarketClearing, Ploadfc, HHourCap, collect(1:1))
	println(DAdecision)
	println(DAcost)
	RTcost = RTNminus1proxy(N, Dn, D, Gn, Gtype, G, Plim, RD, RU, cost, beta, X, fmax, STfmax, L, BranchOutRate, BranchSchOut, BranchForceOut, DAdecision, Pload, HHourCap, VoLL, collect(1:24))
	println(RTcost)

end
