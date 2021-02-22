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
        result["termination_status"] = result_n["termination_status"]
        result["primal_status"] = result_n["primal_status"]
        result["dual_status"] = result_n["dual_status"]

        # TODO: Use objective-based feasibility test.
        if result_n["objective"] > 1.0e-6 || result_n["primal_status"] !== FEASIBLE_POINT
            result["primal_status"] = INFEASIBLE_POINT
        end

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

"Simulates a network dataset, `data`, with controls in `result`, by solving a sequence of convex programs using `optimizer`."
function simulate!(data::Dict{String, <:Any}, schedules, result, result_mn, optimizer)
    data_tmp = deepcopy(data)

    # Ensure data dictionaries have the required properties.
    @assert WM._IM.ismultinetwork(data_tmp) == false
    @assert haskey(data_tmp, "time_series") && haskey(data_tmp["time_series"], "num_steps")

    # Add control time series information from `result` to `data`.
    _update_control_time_series!(data_tmp, schedules, result)

    for n in 1:data["time_series"]["num_steps"]
        # Load and fix data at the current time step.
        WM._IM.load_timepoint!(data_tmp, n)
        WM.fix_all_indicators!(data_tmp)

        # Solve the single time period subproblem.
        result_n = WM.solve_wf(data_tmp, CDWaterModel, optimizer; relax_integrality = true)
        result_mn["termination_status"] = result_n["termination_status"]
        result_mn["primal_status"] = result_n["primal_status"]
        result_mn["dual_status"] = result_n["dual_status"]

        # Update solution data with the recent solution.
        WM._IM.update_data!(result_mn["solution"]["nw"][string(n)], result_n["solution"])

        if result_n["objective"] > 1.0e-6 || result_n["primal_status"] !== FEASIBLE_POINT
            result_mn["primal_status"] = INFEASIBLE_POINT
        end

        # If the problem is not feasible, exit the loop.
        if result_mn["primal_status"] !== FEASIBLE_POINT
            result_mn["last_nw"] = n
            break # Exits the main simulation loop.
        else
            result_mn["last_nw"] = nothing
        end

        # Update tank heads for the next time step using the current solution.
        _update_initial_tank_heads!(data_tmp, result_n)
    end

    if result_mn["primal_status"] === FEASIBLE_POINT
        # Check if the tank levels are recovered at the end of the period.
        nw_ids = sort([parse(Int, x) for x in keys(result_mn["solution"]["nw"])])
        tank_sol_start = result_mn["solution"]["nw"][string(nw_ids[1])]["tank"]
        tank_sol_end = result_mn["solution"]["nw"][string(nw_ids[end])]["tank"]

        @assert Set(keys(tank_sol_start)) == Set(keys(tank_sol_end))
        tank_ids = sort(collect(keys(tank_sol_start)))

        recovered = all(tank_sol_end[i]["V"] >= tank_sol_start[i]["V"] for i in tank_ids)

        if !recovered
            result_mn["last_nw"] = nw_ids[end]
            result_mn["primal_status"] = INFEASIBLE_POINT
        end
    end
end


function _update_control_time_series!(data::Dict{String, <:Any}, schedule, solution)
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

        statuses = [solution[n][k] for n in sort(collect(keys(schedule)))]
        data["time_series"][comp_type][comp_id]["status"] = Array{Int64, 1}(statuses)
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


function _set_median_tank_heads!(data::Dict{String, <:Any})
    for (i, tank) in data["tank"]
        mid_level = tank["min_level"] + 0.5 * (tank["max_level"] - tank["min_level"])
        tank["init_level"] = mid_level
    end
end


function _update_initial_tank_heads!(data::Dict{String, <:Any}, result::Dict{String, <:Any})
    for (i, tank) in data["tank"]
        dV = result["solution"]["tank"][i]["q"] * data["time_series"]["time_step"]
        tank["init_level"] -= dV / (0.25 * pi * tank["diameter"]^2)
    end
end