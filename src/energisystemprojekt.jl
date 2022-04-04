# I den här filen bygger ni modellen. Notera att det är skrivet som en modul, dvs ett paket.
# Så när ni ska använda det, så skriver ni Using energisystemprojekt i er REPL, då får ni ut det ni
# exporterat. Se rad 9.

module energisystemprojekt

using JuMP, AxisArrays, Gurobi, UnPack, Plots

export runmodel

include("input_energisystemprojekt.jl")

function buildmodel(input)

    println("\nBuilding model...")

    @unpack REGION, PLANT, HOUR, numregions, load, maxcap, annualcost, runcost, fucost, hydro_inflow, pv_cf, wind_cf, ef_gas = input

    m = Model(Gurobi.Optimizer)

    @variables m begin

        Electricity[r in REGION, p in PLANT, h in HOUR]       >= 0        # MWh/h
        Capacity[r in REGION, p in PLANT]                     >= 0        # MW
        Reservoar[h in HOUR]                                  >= 0        # MWh
        Systemcost[r in REGION]                               >= 0        # Euro
        Emissions[r in REGION]                                >= 0        # Ton CO2

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
            sum(Electricity[r, p, h] for p in PLANT) >= load[r,h]
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

    end #constraints


    @objective m Min begin
        sum(Systemcost[r] for r in REGION)
    end # objective

    return (;m, Capacity, Emissions, Electricity)

end # buildmodel

function runmodel()

    input = read_input()

    model = buildmodel(input)

    @unpack m, Capacity, Emissions, Electricity = model

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

    println("Cost (M€): ", Cost_result)
    println(Capacity_result)
    println(emissions_result, " ton CO2")

#    for r in REGION
#        y = sum(Capacity[p ,r] for p in PLANT)
#        x = HOUR[1,length(HOUR)]
#        plot(x,y)
#    end
# start of plot for exercise 1

    x = 147:160#651
    y_w = value.(Electricity[:DE, :Wind, x])
    y_pv = zeros(length(x))
    y_gas = zeros(length(x))
    for z in 1:length(x)
        y_pv[z] = value.(Electricity[:DE, :Wind, z]) + value.(Electricity[:DE, :PV, z])
        y_gas[z] = y_pv[z] + value.(Electricity[:DE, :Gas, z])
    end

    plot(x,y_gas)
    plot(x,y_pv)
    plot(x,y_w)
    nothing

end #runmodel


end # module
