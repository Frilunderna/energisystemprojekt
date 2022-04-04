# I den här filen kan ni stoppa all inputdata.
# Läs in datan ni fått som ligger på Canvas genom att använda paketen CSV och DataFrames

using CSV, DataFrames

function read_input()
println("\nReading Input Data...")
folder = dirname(@__FILE__)

#Sets
REGION = [:DE, :SE, :DK]
PLANT = [:Hydro, :Gas, :PV, :Wind] # Add all plants
HOUR = 1:8760

#Parameters
numregions = length(REGION)
numhours = length(HOUR)

timeseries = CSV.read("$folder\\TimeSeries.csv", DataFrame)
wind_cf = AxisArray(ones(numregions, numhours), REGION, HOUR)
load = AxisArray(zeros(numregions, numhours), REGION, HOUR)
pv_cf = AxisArray(ones(numregions, numhours), REGION, HOUR)
hydro_inflow = AxisArray(timeseries[:,"Hydro_inflow"], HOUR)


    for r in REGION
        wind_cf[r, :]=timeseries[:, "Wind_"*"$r"]                             # 0-1, share of installed cap
        pv_cf[r, :]=timeseries[:, "PV_"*"$r"]
        load[r, :]=timeseries[:, "Load_"*"$r"]                                # [MWh]
    end


maxcaptable = [                                                             # GW
        # PLANT      DE             SE              DK
        :Hydro       0              14              0
        :Gas         10000          10000           10000
        :PV          460            75              60
        :Wind        180            280             90
        ]

maxcap = AxisArray(maxcaptable[:,2:end]'.*1000, REGION, PLANT) # MW

investmentcost = [              #euro/KW
        # Plant
        :Hydro          0
        :Gas            550
        :PV             600
        :Wind           1100
        ]
invcost = AxisArray(vec(investmentcost[:,2]'./1000), PLANT) #euro/MW

runningcost = [                 #euro/MWh
        # Plant
        :Hydro          0.1
        :Gas            2
        :PV             0.1
        :Wind           0.1
        ]
runcost = AxisArray(vec(runningcost[:,2]), PLANT) #euro/MWh

fuelcost = [                     #euro/MWh_el
        :Hydro          0
        :Gas            22/0.4
        :PV             0
        :Wind           0
        ]
fucost = AxisArray(vec(fuelcost[:,2]), PLANT) #euro/MWh_el

lifetime = [
        :Hydro          80
        :Gas            30
        :PV             25
        :Wind           25
        ]

lt = AxisArray(vec(lifetime[:,2]), PLANT) #years

efficiency =  [
        :Hydro          1
        :Gas            0.4
        :PV             1
        :Wind           1
        ]
eff = AxisArray(vec(efficiency[:,2]), PLANT) #MWh/MWh

ef_gas = 0.202

discountrate = 0.05

annualcost= AxisArray(ones(length(PLANT)), PLANT)
for p in PLANT
        annualcost[p]=invcost[p]*discountrate/(1-(1/((1+discountrate)^lt[p])))
end


      return (; REGION, PLANT, HOUR, numregions, load, maxcap, annualcost, runcost, fucost, hydro_inflow, pv_cf, wind_cf, ef_gas)

end # read_input
