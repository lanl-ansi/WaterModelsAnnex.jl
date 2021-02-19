"Simulates a network dataset with file path `path`, by solving a sequence of convex programs using `optimizer`."
function simulate(path::String, optimizer)
     result = simulate(WM.parse_file(path), optimizer)
end


"Simulates a network dataset, `data`, with controls in `result`, by solving a sequence of convex programs using `optimizer`."
function simulate!(data::Dict{String, <:Any}, result::Dict{String, <:Any}, optimizer)
    # Ensure data dictionaries have the required properties.
    @assert WM._IM.ismultinetwork(data) == false
    @assert haskey(data, "time_series") && haskey(data["time_series"], "num_steps")

    # Add control time series information from `result` to `data`.
    _update_control_time_series!(data, result)

    for n in 1:data["time_series"]["num_steps"]
        # Load and fix data at the current time step.
        WM._IM.load_timepoint!(data, n)
        WM.fix_all_indicators!(data)

        # Solve the single time period subproblem.
        result_n = WM.solve_wf(data, CDWaterModel, optimizer; relax_integrality = true)

        # Update solution data with the recent solution.
        WM._IM.update_data!(result["solution"]["nw"][string(n)], result_n["solution"])

        # If the problem is not feasible, exit the loop.
        if result_n["primal_status"] !== FEASIBLE_POINT
            result["termination_status"] = result_n["termination_status"]
            result["primal_status"] = result_n["primal_status"]
            result["dual_status"] = result_n["dual_status"]
            break # Exits the main simulation loop.
        end

        # Update tank heads for the next time step using the current solution.
        _update_initial_tank_heads!(data, result_n)
    end

    return result
end


function _update_control_time_series!(data::Dict{String, <:Any}, result::Dict{String, <:Any})
    if !haskey(data["time_series"], "pump")
        data["time_series"]["pump"] = Dict{String, Any}()
    end

    if !haskey(data["time_series"], "valve")
        data["time_series"]["valve"] = Dict{String, Any}()
    end

    nw_ids = sort([parse(Int, x) for x in keys(result["solution"]["nw"])])

    for (i, pump) in data["pump"]
        sol = result["solution"]["nw"]

        if !haskey(data["time_series"]["pump"], i)
            data["time_series"]["pump"][i] = Dict{String, Any}()
        end

        statuses = Array{Int64, 1}([Int(round(sol[string(n)]["pump"][i]["status"])) for n in nw_ids])
        data["time_series"]["pump"][i]["status"] = Array{Int64, 1}(statuses)
    end

    for (i, valve) in data["valve"]
        sol = result["solution"]["nw"]

        if !haskey(data["time_series"]["valve"], i)
            data["time_series"]["valve"][i] = Dict{String, Any}()
        end

        statuses = Array{Int64, 1}([Int(round(sol[string(n)]["valve"][i]["status"])) for n in nw_ids])
        data["time_series"]["valve"][i]["status"] = statuses
    end
end


function _update_initial_tank_heads!(data::Dict{String, <:Any}, result::Dict{String, <:Any})
    for (i, tank) in data["tank"]
        dV = result["solution"]["tank"][i]["q"] * data["time_series"]["time_step"]
        tank["init_level"] -= 4.0 * dV / (pi * tank["diameter"]^2)
    end
end