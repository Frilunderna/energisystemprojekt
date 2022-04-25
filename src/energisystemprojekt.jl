# I den här filen bygger ni modellen. Notera att det är skrivet som en modul, dvs ett paket.
# Så när ni ska använda det, så skriver ni Using energisystemprojekt i er REPL, då får ni ut det ni
# exporterat. Se rad 9.

module energisystemprojekt

using JuMP, AxisArrays, Gurobi, UnPack, Plots, StatsPlots

export runmodel

include("input_energisystemprojekt.jl")

function buildmodel(input)

    println("\nBuilding model...")

    @unpack REGION, PLANT, HOUR, numregions, load, maxcap, annualcost, runcost, fucost, hydro_inflow, pv_cf, wind_cf, ef_gas, eff = input

    m = Model(Gurobi.Optimizer)

    @variables m begin

        Electricity[r in REGION, p in PLANT, h in HOUR]       >= 0        # MWh/h
        Capacity[r in REGION, p in PLANT]                     >= 0        # MW
        Reservoar[h in HOUR]                                  >= 0        # MWh
        Systemcost[r in REGION]                               >= 0        # Euro
        Emissions[r in REGION]                                >= 0        # Ton CO2
        Battery_Storage[r in REGION, h in HOUR]               >= 0        # MWh
        Transmission[r1 in REGION, r2 in REGION, h in HOUR]   >= 0        # MWh
                    #From region   To region
        Capacity_T[r1 in REGION, r2 in REGION]                >= 0        # MW

    end #variables


    #Variable bounds
    for r in REGION, p in PLANT
        set_upper_bound(Capacity[r, p], maxcap[r, p])
    end
    for h in HOUR
        set_upper_bound(Reservoar[h],33*1000000)
    end


    @constraints m begin
        Generation_hydro[r in REGION, h in HOUR],
            Electricity[r, :Hydro, h] <= Capacity[r, :Hydro]
        Generation_gas[r in REGION, h in HOUR],
            Electricity[r, :Gas, h] <= Capacity[r, :Gas]
        Generation_pv[r in REGION, h in HOUR],
            Electricity[r, :PV, h] <= Capacity[r, :PV] * pv_cf[r, h]
        Generation_wind[r in REGION, h in HOUR],
            Electricity[r, :Wind, h] <= Capacity[r, :Wind] * wind_cf[r, h]
        Demand[r in REGION, h in HOUR],
            sum(Electricity[r, p, h] for p in PLANT) >= load[r,h] #- sum(Transmission[r, r2, h] for r2 in REGION) - Battery_Storage[r, h]
        SystemCost[r in REGION],
            Systemcost[r] >= sum(Capacity[r, p]*annualcost[p] for p in PLANT) + sum(Electricity[r, p, h] * (runcost[p] + fucost[p]) for p in PLANT, h in HOUR) # sum of all annualized costs
        Reserv1[h in 2:length(HOUR)],
            Reservoar[h] == Reservoar[h-1]+hydro_inflow[h]-Electricity[:SE,:Hydro,h]
        Reserv2[h in HOUR],
            sum(Electricity[r,:Hydro,h] for r in REGION) <= Reservoar[h]
        Reserv3,
            Reservoar[1] == Reservoar[length(HOUR)]
        Gas_emissions[r in REGION],
            Emissions[r] == sum(Electricity[r, :Gas, h] for h in HOUR) * ef_gas
        #Start of exercise 2
        # a)
        Joint_Emissions,
            sum(Emissions[r] for r in REGION) <= (8552556 + 4978326 + 125243664)*0.1
        # b)
        Batteries[r in REGION, h in HOUR],
            Battery_Storage[r, h] <= Capacity[r, :Battery]
        Generation_battery[r in REGION, h in HOUR],
            Electricity[r, :Battery, h] <= Battery_Storage[r, h] * 0.9
        Battery_change[r in REGION, h in 2:length(HOUR)],
            Battery_Storage[r, h] == Battery_Storage[r, h-1] - Electricity[r, :Battery, h]/0.9 + sum(Electricity[r, p, h] for p in PLANT) - load[r, h] - sum(Transmission[r, r2, h] for r2 in REGION)
        Batteries2[r in REGION],
            Battery_Storage[r,1] == Battery_Storage[r,length(HOUR)]
        # Start of exercise 3
        Transmission_in[h in HOUR, r in REGION],
            Electricity[r, :Transmission, h] <= sum(Transmission[r2, r, h] for r2 in REGION) * 0.98
        Transmission_Capa[r in REGION],
            Capacity[r, :Transmission] == sum(Capacity_T[r2,r] for r2 in REGION)
        Transmission_Cap[r1 in REGION, r2 in REGION, h in HOUR],
            Transmission[r1,r2,h] <= Capacity_T[r1,r2]
        Transmission_equal[r1 in REGION, r2 in REGION],
            Capacity_T[r1,r2] == Capacity_T[r2,r1]
        Transmission_same[r in REGION],
            Capacity_T[r,r] == 0
        # Exercise 4
        Generation_nuclear[r in REGION, h in HOUR],
            Electricity[r, :Nuclear, h] <= Capacity[r, :Nuclear]


    end #constraints


    @objective m Min begin
        sum(Systemcost[r] for r in REGION)
    end # objective

    return (;m, Capacity, Emissions, Electricity, HOUR, load, Capacity_T)

end # buildmodel

function runmodel()

    input = read_input()

    model = buildmodel(input)

    @unpack m, Capacity, Emissions, Electricity, HOUR, load, Capacity_T = model

    println("\nSolving model...")

    status = optimize!(m)


    if termination_status(m) == MOI.OPTIMAL
        println("\nSolve status: Optimal")
    elseif termination_status(m) == MOI.TIME_LIMIT && has_values(m)
        println("\nSolve status: Reached the time-limit")
    else
        error("The model was not solved correctly.")
    end

    Cost_result = objective_value(m)/1000000 # M€
    Capacity_result = value.(Capacity)
    emissions_result = value.(Emissions)
    el=value.(Electricity)
    capacity_t = value.(Capacity_T)

    println("Cost (M€): ", Cost_result)
    println(Capacity_result)
    println(emissions_result, " ton CO2")
    println(capacity_t)


    x = value.(147:651)
    y_w = zeros(length(x))
    y_pv = zeros(length(x))
    y_gas = zeros(length(x))
    y_battery = zeros(length(x))
    y_transmission = zeros(length(x))
    y_nuclear = zeros(length(x))
    y_load = zeros(length(x))
    for z in 1:length(x)
        y_w[z] = value.(Electricity[:DE, :Wind, z+146])
        y_pv[z] = y_w[z] + value.(Electricity[:DE, :PV, z+146])
        y_gas[z] = y_pv[z] + value.(Electricity[:DE, :Gas, z+146])
        y_battery[z] = y_gas[z] + value.(Electricity[:DE, :Battery, z+146])
        y_transmission[z] = y_battery[z] + value.(Electricity[:DE, :Transmission, z+146])
        y_nuclear[z] = y_transmission[z] + value.(Electricity[:DE, :Nuclear, z+146])
        y_load[z] = value.(load[:DE, z+146])
    end


    display(plot(x, [y_nuclear, y_transmission, y_gas, y_pv, y_w, y_load], label=["gas, sol, vind, batteri, transmission och kärnkraft" "gas, sol, vind, batteri och transmission" "gas, sol, vind och batteri" "gas, sol och vind" "sol och vind"  "vind" "load"]))
    savefig("Exercise4_tot.pdf")

    wind = [value.(Capacity[:DE,:Wind]), value.(Capacity[:SE,:Wind]), value.(Capacity[:DK,:Wind])]
    pv = [value.(Capacity[:DE,:PV]),value.(Capacity[:SE,:PV]),value.(Capacity[:DK,:PV])]
    gas = [value.(Capacity[:DE,:Gas]), value.(Capacity[:SE,:Gas]), value.(Capacity[:DK,:Gas])]
    hydro = [value.(Capacity[:DE,:Hydro]), value.(Capacity[:SE,:Hydro]), value.(Capacity[:DK,:Hydro])]
    battery = [value.(Capacity[:DE,:Battery]), value.(Capacity[:SE,:Battery]), value.(Capacity[:DK,:Battery])]
    transmission = [value.(Capacity[:DE,:Transmission]), value.(Capacity[:SE,:Transmission]), value.(Capacity[:DK,:Transmission])]
    nuclear = [value.(Capacity[:DE,:Nuclear]), value.(Capacity[:SE,:Nuclear]), value.(Capacity[:DK,:Nuclear])]

    ticklabel = ["DE" "SE" "DK"]
    display(
    groupedbar([wind pv gas hydro battery transmission nuclear],
    bar_position = :stack,
    bar_width=0.7,
    xticks=(1:3, ticklabel),
    label=["Wind" "PV" "Gas" "Hydro" "Battery" "Transmission" "Nuclear"])
    )
    savefig("Exercise4_totcap.pdf")

    wind = [value.(sum(Electricity[:DE,:Wind, h] for h in HOUR)), value.(sum(Electricity[:SE,:Wind, h] for h in HOUR)), value.(sum(Electricity[:DK,:Wind, h] for h in HOUR))]
    pv = [value.(sum(Electricity[:DE,:PV, h] for h in HOUR)),value.(sum(Electricity[:SE,:PV, h] for h in HOUR)),value.(sum(Electricity[:DK,:PV, h] for h in HOUR))]
    gas = [value.(sum(Electricity[:DE,:Gas, h] for h in HOUR)), value.(sum(Electricity[:SE,:Gas, h] for h in HOUR)), value.(sum(Electricity[:DK,:Gas, h] for h in HOUR))]
    hydro = [value.(sum(Electricity[:DE,:Hydro, h] for h in HOUR)), value.(sum(Electricity[:SE,:Hydro, h] for h in HOUR)), value.(sum(Electricity[:DK,:Hydro, h] for h in HOUR))]
    battery = [value.(sum(Electricity[:DE,:Battery, h] for h in HOUR)), value.(sum(Electricity[:SE,:Battery, h] for h in HOUR)), value.(sum(Electricity[:DK,:Battery, h] for h in HOUR))]
    transmission = [value.(sum(Electricity[:DE,:Transmission, h] for h in HOUR)), value.(sum(Electricity[:SE,:Transmission, h] for h in HOUR)), value.(sum(Electricity[:DK,:Transmission, h] for h in HOUR))]
    nuclear = [value.(sum(Electricity[:DE,:Nuclear, h] for h in HOUR)), value.(sum(Electricity[:SE,:Nuclear, h] for h in HOUR)), value.(sum(Electricity[:DK,:Nuclear, h] for h in HOUR))]

    ticklabel = ["DE" "SE" "DK"]
    display(
    groupedbar([wind pv gas hydro battery transmission nuclear],
    bar_position = :stack,
    bar_width=0.7,
    xticks=(1:3, ticklabel),
    label=["Wind" "PV" "Gas" "Hydro" "Battery" "Transmission" "Nuclear"])
    )
    savefig("Exercise4_totele.pdf")

    nothing

end #runmodel


end # module
