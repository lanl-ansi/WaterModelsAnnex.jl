"Simulates a network dataset with file path `path`, by solving a sequence of convex programs using `optimizer`."
function simulate(path::String, optimizer)
     result = simulate(WM.parse_file(path), optimizer)
end


"Simulates a network dataset, `data`, with controls in `result`, by solving a sequence of convex programs using `optimizer`."
function simulate!(data::Dict{String, <:Any}, result::Dict{String, <:Any}, optimizer)
    data_tmp = deepcopy(data)
    result_tmp = deepcopy(result)

    # Ensure data dictionaries have the required properties.
    @assert WM._IM.ismultinetwork(data_tmp) == false
    @assert haskey(data_tmp, "time_series") && haskey(data_tmp["time_series"], "num_steps")

    # Add control time series information from `result` to `data`.
    _update_control_time_series!(data_tmp, result_tmp)
    println(data_tmp["time_series"]["pump"])

    for n in 1:data_tmp["time_series"]["num_steps"]
        # Load and fix data at the current time step.
        WM._IM.load_timepoint!(data_tmp, n)
        WM.fix_all_indicators!(data_tmp)
        println(n, " ", data_tmp["tank"]["1"]["init_level"], " ", data_tmp["tank"]["1"]["max_level"])

        # Solve the single time period subproblem.
        result_n = WM.solve_wf(data_tmp, CDWaterModel, optimizer; relax_integrality = true)
        result_tmp["termination_status"] = result_n["termination_status"]
        result_tmp["primal_status"] = result_n["primal_status"]
        result_tmp["dual_status"] = result_n["dual_status"]

        println(result_n["termination_status"], " ", result_n["objective"], " ", result_n["primal_status"])

        # TODO: Use objective-based feasibility test.
        if result_n["objective"] > 1.0e-6 || !(result_n["primal_status"] in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT])
            result_tmp["primal_status"] = INFEASIBLE_POINT
        end

        # Update solution data with the recent solution.
        WM._IM.update_data!(result_tmp["solution"]["nw"][string(n)], result_n["solution"])

        # If the problem is not feasible, exit the loop.
        if !(result_n["primal_status"] in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT])
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

        sol_nw = result_mn["solution"]["nw"][string(n)]

        if haskey(sol_nw, "pump")
            for (a, pump) in sol_nw["pump"]
                if haskey(result_n["solution"], "pump") && haskey(result_n["solution"]["pump"], a)
                    continue
                elseif haskey(result_n["solution"], "pump") && !haskey(result_n["solution"]["pump"], a)
                    pump["status"] = 0
                    # for (key, value) in pump
                    #     pump[key] = 0.0
                    # end
                else
                    pump["status"] = 0
                    # for (key, value) in pump
                    #     pump[key] = 0.0
                    # end
                end
            end
        end

        if haskey(sol_nw, "valve")
            for (a, valve) in sol_nw["valve"]
                if haskey(result_n["solution"], "valve") && haskey(result_n["solution"]["valve"], a)
                    continue
                elseif haskey(result_n["solution"], "valve") && !haskey(result_n["solution"]["valve"], a)
                    valve["status"] = 0
                    # for (key, value) in valve
                    #     valve[key] = 0.0
                    # end
                else
                    valve["status"] = 0
                    # for (key, value) in valve
                    #     valve[key] = 0.0
                    # end
                end
            end
        end

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

        if n < data["time_series"]["num_steps"]
            # Update tank heads for the next time step using the current solution.
            _update_initial_tank_heads!(data_tmp, result_n)
        end
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

        statuses = Array{Int64, 1}([Int(round(sol[string(n)]["pump"][i]["status"])) for n in nw_ids])
        data["time_series"]["pump"][i]["status"] = Array{Int64, 1}(statuses)
        data["time_series"]["pump"][i]["z_min"] = Array{Int64, 1}(statuses)
        data["time_series"]["pump"][i]["z_max"] = Array{Int64, 1}(statuses)
    end

    for (i, valve) in data["valve"]
        sol = result["solution"]["nw"]

        if !haskey(data["time_series"]["valve"], i)
            data["time_series"]["valve"][i] = Dict{String, Any}()
        end

        statuses = Array{Int64, 1}([Int(round(sol[string(n)]["valve"][i]["status"])) for n in nw_ids])
        data["time_series"]["valve"][i]["status"] = statuses
        data["time_series"]["valve"][i]["z_min"] = statuses
        data["time_series"]["valve"][i]["z_max"] = statuses
    end
end


function _update_tank_time_series!(data::Dict{String, <:Any}, result::Dict{String, <:Any})
    if !haskey(data["time_series"], "tank")
        data["time_series"]["tank"] = Dict{String, Any}()
    end

    nw_ids = sort([parse(Int, x) for x in keys(result["solution"]["nw"])])

    for (i, tank) in data["tank"]
        sol = result["solution"]["nw"]

        if !haskey(data["time_series"]["tank"], i)
            data["time_series"]["tank"][i] = Dict{String, Any}()
        end

        V = Array{Int64, 1}([Int(round(sol[string(n)]["tank"][i]["V"])) for n in nw_ids])
        surface_area = 0.25 * pi * tank["diameter"]^2
        init_level = V ./ surface_area

        data["time_series"]["tank"][i]["init_level"] = init_level
    end
end


function _set_median_tank_heads!(data::Dict{String, <:Any})
    for (i, tank) in data["tank"]
        level_delta = 0.5 * (tank["max_level"] - tank["min_level"])
        tank["init_level"] = tank["min_level"] + level_delta
    end
end


function _update_initial_tank_heads!(data::Dict{String, <:Any}, result::Dict{String, <:Any})
    for (i, tank) in data["tank"]
        dV = result["solution"]["tank"][i]["q"] * data["time_series"]["time_step"]
        tank["init_level"] -= dV / (0.25 * pi * tank["diameter"]^2)
    end
end