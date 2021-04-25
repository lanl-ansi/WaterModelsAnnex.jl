"Simulates a network dataset with file path `path`, by solving a sequence of convex programs using `optimizer`."
function simulate(path::String, optimizer)
     result = simulate(WM.parse_file(path), optimizer)
end


"Simulates a network dataset `data`, by solving a sequence of convex programs using `optimizer`."
function simulate!(data::Dict{String, <:Any}, optimizer::Any)
    for time_step_index in 1:data["time_series"]["num_steps"]
        WM._IM.load_timepoint!(data, time_step_index)
        result = WM.solve_wf(data, CDWaterModel, optimizer)
        update_tank_heads!(data, result["solution"])
    end
end


"Simulates a network dataset, `data`, with controls in `result`, by solving a sequence of convex programs using `optimizer`."
function simulate!(data::Dict{String, <:Any}, result::Dict{String, <:Any}, optimizer)
    data_tmp = deepcopy(data)
    result_tmp = deepcopy(result)

    # Ensure data dictionaries have the required properties.
    @assert !WM._IM.ismultinetwork(data_tmp)
    @assert haskey(data_tmp, "time_series")

    # Add control time series information from `result` to `data`.
    _update_control_time_series!(data_tmp, result_tmp)

    for n in 1:data_tmp["time_series"]["num_steps"]
        # Load and fix data at the current time step.
        WM._IM.load_timepoint!(data_tmp, n)
        WM.fix_all_indicators!(data_tmp)

        # Solve the single time period subproblem.
        result_n = WM.solve_wf(data_tmp, CDWaterModel, optimizer; relax_integrality = true)
        result_tmp["termination_status"] = result_n["termination_status"]
        result_tmp["primal_status"] = result_n["primal_status"]
        result_tmp["dual_status"] = result_n["dual_status"]

        # TODO: Use objective-based feasibility test.
        if result_n["objective"] > 1.0e-6 || !(result_n["primal_status"] in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT])
            result_tmp["primal_status"] = INFEASIBLE_POINT
        end

        # Update solution data with the recent solution.
        WM._IM.update_data!(result_tmp["solution"]["nw"][string(n)], result_n["solution"])

        # If the problem is not feasible, exit the loop.
        if !feasible_simulation_result(result_n)
            result_tmp["termination_status"] = result_n["termination_status"]
            result_tmp["primal_status"] = result_n["primal_status"]
            result_tmp["dual_status"] = result_n["dual_status"]
            break # Exits the main simulation loop.
        end

        # Update tank heads for the next time step using the current solution.
        _update_initial_tank_heads!(data_tmp, result_n)
    end

    return result_tmp
end

"Simulates a network dataset, `data`, with controls in `result`, by solving a sequence of convex programs using `optimizer`."
function simulate!(data::Dict{String, <:Any}, schedules, result, result_mn, optimizer)
    data_tmp = deepcopy(data)

    # Ensure data dictionaries have the required properties.
    @assert WM._IM.ismultinetwork(data_tmp) == false
    @assert haskey(data_tmp, "time_series")
    @assert haskey(data_tmp["time_series"], "num_steps")

    # Add control time series information from `result` to `data`.
    _update_control_time_series!(data_tmp, schedules, result)

    for n in 1:data["time_series"]["num_steps"] - 1
        # Load and fix data at the current time step.
        WM._IM.load_timepoint!(data_tmp, n)
        WM.fix_all_indicators!(data_tmp)

         # Solve the single time period subproblem.
        result_n = WM.solve_wf(data_tmp, CDWaterModel, optimizer; relax_integrality = true)
        result_mn["termination_status"] = result_n["termination_status"]
        result_mn["primal_status"] = result_n["primal_status"]
        result_mn["dual_status"] = result_n["dual_status"]
        sol_nw = result_mn["solution"]["nw"][string(n)]
        WM._IM.update_data!(sol_nw, result_n["solution"])

        if haskey(sol_nw, "pump")
            for (a, pump) in sol_nw["pump"]
                if haskey(result_n["solution"], "pump") && haskey(result_n["solution"]["pump"], a)
                    pump["status"] = WM.STATUS(Int(round(result_n["solution"]["pump"][a]["status"])))
                elseif haskey(result_n["solution"], "pump") && !haskey(result_n["solution"]["pump"], a)
                    pump["status"] = WM.STATUS_INACTIVE
                else
                    pump["status"] = WM.STATUS_INACTIVE
                end
            end
        end

        if haskey(sol_nw, "valve")
            for (a, valve) in sol_nw["valve"]
                if haskey(result_n["solution"], "valve") && haskey(result_n["solution"]["valve"], a)
                    valve["status"] = WM.STATUS(Int(round(result_n["solution"]["valve"]["status"])))
                elseif haskey(result_n["solution"], "valve") && !haskey(result_n["solution"]["valve"], a)
                    valve["status"] = WM.STATUS_INACTIVE
                else
                    valve["status"] = WM.STATUS_INACTIVE
                end
            end
        end

        if !feasible_simulation_result(result_n)
            result_mn["primal_status"] = INFEASIBLE_POINT
            result_mn["objective"] = Inf
        else
            result_mn["primal_status"] = FEASIBLE_POINT
            result_mn["objective"] = 0.0
        end

        # If the problem is not feasible, exit the loop.
        if !feasible_simulation_result(result_mn)
            result_mn["primal_status"] = INFEASIBLE_POINT
            result_mn["objective"] = Inf
            result_mn["last_nw"] = n
            break # Exits the main simulation loop.
        else
            result_mn["primal_status"] = FEASIBLE_POINT
            result_mn["objective"] = 0.0
            result_mn["last_nw"] = nothing
        end

        # Update tank heads for the next time step using the current solution.
        _update_initial_tank_heads!(data_tmp, result_n, n+1)
    end

    if feasible_simulation_result(result_mn)
        _update_final_tank_heads!(data_tmp, result_mn)
    end

    if feasible_simulation_result(result_mn)
        # Check if the tank levels are recovered at the end of the period.
        node_indices = [string(x["node"]) for (i, x) in data["tank"]]
        nw_ids = sort([parse(Int, x) for x in keys(result_mn["solution"]["nw"])])
        node_sol_start = result_mn["solution"]["nw"][string(nw_ids[1])]["node"]
        node_sol_end = result_mn["solution"]["nw"][string(nw_ids[end])]["node"]

        @assert Set(keys(node_sol_start)) == Set(keys(node_sol_end))
        recovered = all([node_sol_end[i]["h"] >= node_sol_start[i]["h"] for i in node_indices])

        if !recovered
            result_mn["last_nw"] = nw_ids[end-1]
            result_mn["primal_status"] = INFEASIBLE_POINT
            result_mn["objective"] = Inf
        else
            result_mn["primal_status"] = FEASIBLE_POINT
            result_mn["objective"] = 0.0
        end
    end
end


function _update_control_time_series!(data::Dict{String, <:Any}, schedule, result)
    if !haskey(data["time_series"], "pump")
        data["time_series"]["pump"] = Dict{String, Any}()
    end

    if !haskey(data["time_series"], "valve")
        data["time_series"]["valve"] = Dict{String, Any}()
    end

    for (k, var_id) in enumerate(schedule[1][1])
        comp_type = string(var_id.component_type)
        comp_id = string(var_id.component_index)
        comp = data[comp_type][comp_id]

        if !haskey(data["time_series"][comp_type], comp_id)
            data["time_series"][comp_type][comp_id] = Dict{String, Any}()
        end

        statuses = [Int(result[n][k]) for n in sort(collect(keys(schedule)))]
        data["time_series"][comp_type][comp_id]["status"] = Array{Int64, 1}(statuses)
        data["time_series"][comp_type][comp_id]["z_min"] = Array{Int64, 1}(statuses)
        data["time_series"][comp_type][comp_id]["z_max"] = Array{Int64, 1}(statuses)
    end
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

        statuses = Array{Int64, 1}([Int(sol[string(n)]["pump"][i]["status"]) for n in nw_ids])
        data["time_series"]["pump"][i]["status"] = Array{Int64, 1}(statuses)
        data["time_series"]["pump"][i]["z_min"] = Array{Int64, 1}(statuses)
        data["time_series"]["pump"][i]["z_max"] = Array{Int64, 1}(statuses)
    end

    for (i, valve) in data["valve"]
        sol = result["solution"]["nw"]

        if !haskey(data["time_series"]["valve"], i)
            data["time_series"]["valve"][i] = Dict{String, Any}()
        end

        statuses = Array{Int64, 1}([Int(sol[string(n)]["valve"][i]["status"]) for n in nw_ids])
        data["time_series"]["valve"][i]["status"] = statuses
        data["time_series"]["valve"][i]["z_min"] = statuses
        data["time_series"]["valve"][i]["z_max"] = statuses
    end
end


function _update_tank_time_series!(data::Dict{String, <:Any}, result_mn::Dict{String, <:Any})
    if !haskey(data["time_series"], "tank")
        # Initialize the tank time series dictionary.
        data["time_series"]["tank"] = Dict{String, Any}()
    end

    # Get the network indices in the result dictionary.
    nw_ids = sort([parse(Int64, n) for n in keys(result_mn["solution"]["nw"])])

    # Get the relevant piece of the result dictionary.
    sol = result_mn["solution"]["nw"]

    for (i, tank) in data["tank"]
        # Initialize the tank time series dictionary.
        if !haskey(data["time_series"]["tank"], i)
            data["time_series"]["tank"][i] = Dict{String, Any}()
        end

        elev = [data["node"][string(tank["node"])]["elevation"] for n in nw_ids]
        head = [sol[string(n)]["node"][string(tank["node"])]["h"] for n in nw_ids]

        # Store the tank's initial level data in the time series.
        data["time_series"]["tank"][i]["init_level"] = head .- elev
    end
end


function _reset_controllable_component_statuses(data::Dict{String, <:Any})
    for comp_type in ["pump", "regulator", "valve"]
        map(x -> x["status"] = WM.STATUS_UNKNOWN, values(data[comp_type]))
        map(x -> x["z_min"] = 0.0, values(data[comp_type]))
        map(x -> x["z_max"] = 1.0, values(data[comp_type]))
    end
end


function _set_median_tank_heads!(data::Dict{String, <:Any})
    for (i, tank) in data["tank"]
        tank["init_level"] = 0.5 * (tank["max_level"] + tank["min_level"])
    end
end


function _set_median_tank_time_series!(data::Dict{String, <:Any})
    if !haskey(data["time_series"], "tank")
        # Initialize the tank time series dictionary.
        data["time_series"]["tank"] = Dict{String, Any}()
    end

    for (i, tank) in data["tank"]
        # Initialize the tank time series dictionary.
        if !haskey(data["time_series"]["tank"], i)
            data["time_series"]["tank"][i] = Dict{String, Any}()
        end

        # Store the tank's initial level data in the time series.
        mid_level = 0.5 * (tank["min_level"] + tank["max_level"])
        level_array = ones(data["time_series"]["num_steps"], 1) * mid_level
        data["time_series"]["tank"][i]["init_level"] = level_array
    end
end


function update_tank_heads!(data::Dict{String, <:Any}, solution::Dict{String, <:Any})
    for (i, tank) in data["tank"]
        volume_change = solution["tank"][i]["q"] * data["time_series"]["time_step"]
        tank["init_level"] -= volume_change / (0.25 * pi * tank["diameter"]^2)
    end
end


function _update_initial_tank_heads!(data::Dict{String, <:Any}, result::Dict{String, <:Any}, time_index::Int)
    for (i, tank) in data["tank"]
        dV = result["solution"]["tank"][i]["q"] * data["time_series"]["time_step"]
        tank["init_level"] -= dV / (0.25 * pi * tank["diameter"]^2)
        data["time_series"]["tank"][string(i)]["init_level"][time_index] = tank["init_level"]
    end
end


function _update_final_tank_heads!(data::Dict{String, <:Any}, result::Dict{String, <:Any})
    nw_1 = sort([parse(Int, x) for x in collect(keys(result["solution"]["nw"]))])[end-1]
    result_1 = result["solution"]["nw"][string(nw_1)]

    nw_2 = sort([parse(Int, x) for x in collect(keys(result["solution"]["nw"]))])[end]
    result_2 = result["solution"]["nw"][string(nw_2)]

    for (i, tank) in data["tank"]
        dV = result_1["tank"][i]["q"] * data["time_series"]["time_step"]
        #tank["init_level"] -= dV / (0.25 * pi * tank["diameter"]^2)
        #data["time_series"]["tank"][string(i)]["init_level"][nw_2] = tank["init_level"]

        #result_2["tank"][i]["V"] = result_1["tank"][i]["V"] + dV
        result_2["node"][string(tank["node"])]["h"] = 
            result_1["node"][string(tank["node"])]["h"] -
            dV / (0.25 * pi * tank["diameter"]^2)
    end
end


function _update_initial_tank_heads!(data::Dict{String, <:Any}, result::Dict{String, <:Any}, time_step::Float64)
    for (i, tank) in data["tank"]
        dV = result["solution"]["tank"][i]["q"] * time_step
        tank["init_level"] -= dV / (0.25 * pi * tank["diameter"]^2)
    end
end