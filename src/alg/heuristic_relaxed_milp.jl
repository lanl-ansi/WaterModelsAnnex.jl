function build_nw_id_partition(wm::WM.AbstractWaterModel)
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    return [[x] for x in Random.shuffle(nw_ids)]
end


function get_control_settings_at_nw(result::Dict{String, <:Any}, nw::Int)
    sol_nw = result["solution"]["nw"][string(nw)]

    if haskey(sol_nw, "pump")
        pump_ids = sort([parse(Int, x) for x in collect(keys(sol_nw["pump"]))])
        pump_vals = [sol_nw["pump"][string(i)]["status"] for i in pump_ids]
        pump_vars = WM._VariableIndex.(nw, :pump, :z_pump, pump_ids)
    else
        pump_vals, pump_vars = [], []
    end

    if haskey(sol_nw, "valve")
        valve_ids = sort([parse(Int, x) for x in collect(keys(sol_nw["valve"]))])
        valve_vals = [sol_nw["valve"][string(i)]["status"] for i in valve_ids]
        valve_vars = WM._VariableIndex.(nw, :valve, :z_valve, valve_ids)
    else
        valve_vals, valve_vars = [], []
    end

    if haskey(sol_nw, "regulator")
        regulator_ids = sort([parse(Int, x) for x in collect(keys(sol_nw["regulator"]))])
        regulator_vals = [sol_nw["regulator"][string(i)]["status"] for i in regulator_ids]
        regulator_vars = WM._VariableIndex.(nw, :regulator, :z_regulator, regulator_ids)
    else
        regulator_vals, regulator_vars = [], []
    end

    vars = vcat(pump_vars, regulator_vars, valve_vars)
    vals = vcat(pump_vals, regulator_vals, valve_vals)
    return ControlSetting(nw, vars, round.(vals))
end


function get_control_settings_from_result(result::Dict{String, <:Any})
    nw_ids = sort([parse(Int, x) for x in keys(result["solution"]["nw"])])[1:end-1]
    return get_control_settings_at_nw.(Ref(result), nw_ids)
end


function feasibility_pump(wm::WM.AbstractWaterModel)
    objective = JuMP.AffExpr(0.0)

    for nw in sort(collect(WM.nw_ids(wm)))[1:end-1]
        z_pump = vcat(WM.var(wm, nw, :z_pump)...)
        z_regulator = vcat(WM.var(wm, nw, :z_regulator)...)
        z_valve = vcat(WM.var(wm, nw, :z_valve)...)

        z = vcat(z_pump, z_regulator, z_valve)
        z_one = filter(x -> JuMP.is_fixed(x) && round(JuMP.fix_value(x)) == 1.0, z)
        z_zero = filter(x -> JuMP.is_fixed(x) && round(JuMP.fix_value(x)) == 0.0, z)

        if length(z_one) > 0
            objective += (length(z_one) - sum(z_one))
        end

        if length(z_zero) > 0
            objective += sum(z_zero)
        end

        map(x -> JuMP.is_fixed(x) && JuMP.unfix(x), z)
        map(x -> !JuMP.is_binary(x) && JuMP.set_binary(x), z)
     end

    #  JuMP.@constraint(wm.model, objective <= 3.0)
     JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective)
     return WM.optimize_model!(wm)
end


function feasibility_pump_at_nw(wm::WM.AbstractWaterModel, nws_fixed::Vector{Int64}, nws_to_fix::Vector{Int64})
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]

    nw_ids_relaxed = sort(collect(setdiff(Set(nw_ids), Set(vcat(nws_to_fix, nws_fixed)))))
    vars_to_relax = vcat(WM.get_all_binary_vars_at_nw!.(Ref(wm), nw_ids_relaxed)...)
    vars_to_relax_not_y = filter(v -> !occursin("_y", JuMP.name(v)), vars_to_relax)
    map(x -> !JuMP.is_fixed(x) && JuMP.unset_binary(x), vars_to_relax_not_y)

    vars_discrete = vcat(WM.get_all_binary_vars_at_nw!.(
        Ref(wm), vcat(nws_to_fix, nws_fixed))...)
    map(x -> JuMP.set_binary(x), vars_discrete)

    JuMP.optimize!(wm.model) # Solve the relaxed model.

    if JuMP.primal_status(wm.model) === WM._MOI.FEASIBLE_POINT
        vars_to_fix_nx = filter(v -> !occursin("_x", JuMP.name(v)), vars_discrete)
        # vars_to_fix_nxy = filter(v -> !occursin("_y", JuMP.name(v)), vars_to_fix_nx)
        JuMP.fix.(vars_to_fix_nx, round.(JuMP.value.(vars_to_fix_nx)))
        map(x -> JuMP.set_binary(x), vars_to_relax)
        nws_fixed = collect(Set(append!(nws_fixed, nws_to_fix)))
        return true, nws_fixed # Found a feasible solution.
    else
        num_to_unfix = min(length(nws_fixed), 2)
        nw_random_unfix = Random.shuffle(nws_fixed)[1:num_to_unfix]
        vars_to_unfix = vcat(WM.get_all_binary_vars_at_nw!.(Ref(wm), nw_random_unfix)...)
        nws_fixed = collect(setdiff(Set(nws_fixed), Set(nw_random_unfix)))

        map(x -> JuMP.is_fixed(x) && JuMP.unfix(x), vars_to_unfix)
        map(x -> JuMP.set_binary(x), vars_to_unfix)
        map(x -> JuMP.set_binary(x), vars_to_relax)

        return false, nws_fixed
    end
end


function run_feasibility_pump(network_mn::Dict{String, Any}, mip_solver)
    # Specify model options and construct the multinetwork OWF model.
    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, mip_solver)
    feasible, num_iterations = false, 0

    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    nws_to_fix = sort(collect(WM.nw_ids(wm)))[1:end-1]
    nws_fixed = Vector{Int}([])

    while length(nws_fixed) < length(nw_ids) && num_iterations <= 50
        nws_to_fix = collect(Set(setdiff(nw_ids, nws_fixed)))
        num_left = min(length(nws_to_fix), 1)
        nw_to_fix = collect(Random.shuffle(collect(nws_to_fix))[1:num_left])
        feasible, nws_fixed = feasibility_pump_at_nw(wm, nws_fixed, nw_to_fix)
        num_iterations += 1
    end

    binary_vars = filter(v -> JuMP.is_binary(v), JuMP.all_variables(wm.model))
    vars_to_unfix = filter(v -> occursin("_x", JuMP.name(v)), binary_vars)
    vars_to_unfix = filter(v -> occursin("_y", JuMP.name(v)), vars_to_unfix)
    map(v -> JuMP.is_fixed(v) && JuMP.unfix(v), vars_to_unfix)
    map(v -> !JuMP.is_binary(v) && JuMP.set_binary(v), vars_to_unfix)
    result = WM.optimize_model!(wm)

    # if result["primal_status"] !== WM._MOI.FEASIBLE_POINT
    #     result = feasibility_pump(wm)
    # end

    feasible = result["primal_status"] === WM._MOI.FEASIBLE_POINT
    return feasible ? result : nothing
end


function fix_control_setting!(wm::WM.AbstractWaterModel, control_setting::ControlSetting, simulation_result::SimulationResult)
    nw = control_setting.network_id

    for (i, vid) in enumerate(control_setting.variable_indices)
        var = WM._get_variable_from_index(wm, vid)
        !JuMP.is_binary(var) && JuMP.set_binary(var)
        JuMP.is_fixed(var) && JuMP.unfix(var)
        JuMP.fix(var, control_setting.vals[i]; force = true)
    end

    for i in WM.ids(wm, nw, :tank)
        var = WM.var(wm, nw, :q_tank, i)
        JuMP.is_fixed(var) && JuMP.unfix(var)
        JuMP.fix(var, simulation_result.q_tank[i]; force = true)
    end
end


function unfix_control_setting!(wm::WM.AbstractWaterModel, control_setting::ControlSetting)
    nw = control_setting.network_id

    for vid in control_setting.variable_indices
        var = WM._get_variable_from_index(wm, vid)
        JuMP.is_fixed(var) && JuMP.unfix(var)
        !JuMP.is_binary(var) && JuMP.set_binary(var)
    end

    for (i, tank) in WM.ref(wm, nw, :tank)
        var = WM.var(wm, nw, :q_tank, i)
        JuMP.is_fixed(var) && JuMP.unfix(var)
        JuMP.set_lower_bound(var, tank["flow_min"])
        JuMP.set_upper_bound(var, tank["flow_max"])
    end
end


function add_feasibility_cut_from_control_settings!(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting})
    vids = vcat(collect([x.variable_indices for x in control_settings])...)
    vals = vcat(collect([x.vals for x in control_settings])...)
    vars = WM._get_variable_from_index.(Ref(wm), vids)

    one_vars = vars[findall(x -> abs.(round(x)) == 1.0, vals)]
    zero_vars = vars[findall(x -> abs.(round(x)) == 0.0, vals)]
    JuMP.@constraint(wm.model, sum(zero_vars) -
        sum(one_vars) >= 1.0 - length(one_vars))
end


function objective_mean_tank_levels(wm::WM.AbstractWaterModel, nws::Vector{Int})
    objective = JuMP.AffExpr(0.0)

    for nw in nws
        for tank in values(WM.ref(wm, nw, :tank))
            node = WM.ref(wm, nw, :node, tank["node"])
            h = WM.var(wm, nw, :h, tank["node"])
            head_mid = 0.5 * (node["head_min"] + node["head_max"])
            objective += (h - head_mid)^2
        end
    end

    JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective)
end


function objective_predicted_tank_levels(wm::WM.AbstractWaterModel, simulation_results::Vector{SimulationResult}, iteration_num::Int)
    objective = JuMP.AffExpr(0.0)
    tank_heads = Dict{Int, Vector{Float64}}()

    for nw in sort(collect(WM.nw_ids(wm)))[1:end-1]
        time_step = WM.ref(wm, nw, :time_step)

        for tank in values(WM.ref(wm, nw, :tank))
            node = WM.ref(wm, nw, :node, tank["node"])
            q_tank = simulation_results[nw].q_tank[tank["index"]]
            surface_area = 0.25 * pi * tank["diameter"]^2
            h = WM.var(wm, nw + 1, :h, tank["node"])

            if nw == 1
                tank_heads[tank["index"]] = [node["elevation"] + tank["init_level"]]
            end

            push!(tank_heads[tank["index"]],
                tank_heads[tank["index"]][end] -
                q_tank / surface_area * time_step)
        end
    end

    for (i, tank_ts) in tank_heads
        default(show = true)
        plot(1:length(tank_ts), tank_ts, legend=false)

        for nw in sort(collect(WM.nw_ids(wm)))
            tank = WM.ref(wm, nw, :tank, i)
            node = WM.ref(wm, nw, :node, tank["node"])
            min_head = node["elevation"] + tank["min_level"]
            max_head = node["elevation"] + tank["max_level"]

            head_change = 0.0

            if tank_ts[nw] < min_head
                head_change += min_head - tank_ts[nw]
            elseif tank_ts[nw] > max_head
                head_change -= tank_ts[nw] - max_head
            end

            tank_ts[nw] += head_change
            tank_ts[2:end-1] .-= head_change / (length(tank_ts) - 3.0)
        end

        # head_change = (tank_ts[1] - tank_ts[end]) * float(iteration_num)
       
        tank = WM.ref(wm, 1, :tank, i)
        node = WM.ref(wm, 1, :node, tank["node"])
        min_head = node["elevation"] + tank["min_level"]
        max_head = node["elevation"] + tank["max_level"]
        # plot!(1:length(tank_ts), tank_ts, legend=false)
        plot!(1:length(tank_ts), ones(length(tank_ts)) * min_head, legend=false)
        plot!(1:length(tank_ts), ones(length(tank_ts)) * max_head, legend=false)
        readline()
    end

    for nw in sort(collect(WM.nw_ids(wm)))
        for tank in values(WM.ref(wm, nw, :tank))
            h = WM.var(wm, nw, :h, tank["node"])
            objective += (h - tank_heads[tank["index"]][nw])^2
        end
    end

    JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective)
end



function objective_min_losses(wm::WM.AbstractWaterModel, nws::Vector{Int})
    objective = JuMP.AffExpr(0.0)

    for nw in nws
        for pipe in values(WM.ref(wm, nw, :pipe))
            dhp = WM.var(wm, nw, :dhp_pipe, pipe["index"])
            dhn = WM.var(wm, nw, :dhn_pipe, pipe["index"])
            objective += dhp + dhn
        end

        # for pump in values(WM.ref(wm, nw, :pump))
            # objective -= WM.var(wm, nw, :g_pump, pump["index"])
        # end

        # for tank in values(WM.ref(wm, nw, :tank))
        #     node = WM.ref(wm, nw, :node, tank["node"])
        #     h = WM.var(wm, nw, :h, tank["node"])            
        #     head_mid = node["head_min"] + 0.5 * (node["head_max"] - node["head_min"])
        #     objective += (h - head_mid)^2
        # end
    end

    JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective)
end


function constraint_tank_final_range(wm::WM.AbstractWaterModel)
    for nw in sort(collect(WM.nw_ids(wm)))[end]
        for tank in values(WM.ref(wm, nw, :tank))
            node = WM.ref(wm, nw, :node, tank["node"])
            head_min = node["elevation"] + WM.ref(wm, 1, :tank, tank["index"], "init_level")
            head_max = node["elevation"] + WM.ref(wm, nw, :tank, tank["index"], "max_level")
            h = WM.var(wm, nw, :h, tank["node"])

            head_lower = head_min + 0.1 * (head_max - head_min)
            head_upper = head_min + 0.9 * (head_max - head_min)

            JuMP.@constraint(wm.model, h >= head_lower)
            JuMP.@constraint(wm.model, h <= head_upper)
        end
    end
end


function solve_pseudo_relaxation_nw(wm::WM.AbstractWaterModel, nw_random::Int)
    nw_ids = sort(collect(WM.nw_ids(wm)))
    nw_last_binary = findfirst(x -> x == nw_random, nw_ids)
    nw_relaxed = nw_ids[nw_last_binary+1:end]

    vars_to_relax = vcat(WM.get_all_binary_vars_at_nw!.(Ref(wm), nw_relaxed)...)
    map(x -> JuMP.is_fixed(x) && JuMP.unfix(x), vars_to_relax)
    map(x -> !JuMP.is_fixed(x) && JuMP.unset_binary(x), vars_to_relax)
    map(x -> !JuMP.is_fixed(x) && JuMP.set_lower_bound(x, 0.0), vars_to_relax)
    map(x -> !JuMP.is_fixed(x) && JuMP.set_upper_bound(x, 1.0), vars_to_relax)

    result = WM.optimize_model!(wm; relax_integrality = false)
    map(x -> !JuMP.is_binary(x) && JuMP.set_binary(x), vars_to_relax)

    return result
end


function feasibility_pump_set_objective!(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting})
    objective = JuMP.AffExpr(0.0)

    for control_setting in control_settings
        vars = WM._get_variable_from_index.(Ref(wm), control_setting.variable_indices)
        z_one_ids = findall(i -> abs(round(control_setting.vals[i])) == 1.0, 1:length(vars))
        z_zero_ids = findall(i -> abs(round(control_setting.vals[i])) == 0.0, 1:length(vars))
        z_one_vars, z_zero_vars = vars[z_one_ids], vars[z_zero_ids]

        if length(z_one_vars) > 0
            objective += length(z_one_vars) - sum(z_one_vars)
        end

        if length(z_zero_vars) > 0
            objective += sum(z_zero_vars)
        end
    end

    JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective)
end


function feasibility_pump_at_last(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting})
    feasibility_pump_set_objective!(wm, control_settings)
    unfix_control_setting!.(Ref(wm), control_settings)
    return WM.optimize_model!(wm; relax_integrality = false)
end


function run_sequential_heuristic(
    wm::WM.AbstractWaterModel, wm_cq::CQWaterModel,
    result_micp::Dict{String, <:Any}, time_limit::Float64)
    control_settings = get_control_settings_from_result(result_micp)
    map(x -> x.vals = abs.(round.(x.vals)), control_settings)
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    time_elapsed, num_look_back = 0.0, 2

    while time_elapsed < time_limit
        # Simulate the current control schedule forward in time.
        time_elapsed += @elapsed simulation_results =
            simulate_control_settings_sequential(wm_cq, control_settings)

        # If the control schedule is feasible...
        if all(x -> x.feasible, simulation_results)
            # Return the control settings schedule.
            return control_settings
        else # If the control schedule is not feasible...
            # Get the set of times before the first infeasibility.
            nws_feasible = nw_ids[1:max(1, length(simulation_results)-1)]
            WM.Memento.info(LOGGER, "[HEUR] Feasible time indices: $(nws_feasible).")

            # Add a feasibility cut to the WaterModels model up to the infeasible time.
            nws_cut = nw_ids[1:length(simulation_results)]
            add_feasibility_cut_from_control_settings!(wm, control_settings[nws_cut])

            # Pick a random time index at which the controls will be changed.
            num_back = min(num_look_back, length(nws_cut))
            first_id = max(1, length(nws_cut) - num_back)
            nw_random = Random.rand(nws_cut[first_id:end])

            # Unfix the controls at this time step and beyond.
            nws_unfixed = nw_ids[nw_random:end]
            unfix_control_setting!.(Ref(wm), control_settings[nws_unfixed])

            # Fix the controls prior to this time step.
            nws_fixed = nw_random == 1 ? [] : nw_ids[1:max(1, nw_random-1)]
            fix_control_setting!.(Ref(wm), control_settings[nws_fixed], simulation_results[nws_fixed])

            # Solve the corresponding pseudo-relaxation at the random time.
            # objective_min_losses(wm, nws_unfixed[1:end-1])
            objective_mean_tank_levels(wm, sort(collect(WM.nw_ids(wm))))
            time_elapsed += @elapsed result = solve_pseudo_relaxation_nw(wm, nw_random)

            if result["primal_status"] === WM._MOI.FEASIBLE_POINT
                control_settings = get_control_settings_from_result(result)
                map(x -> x.vals = abs.(round.(x.vals)), control_settings)
                num_look_back = 2
            elseif length(nws_unfixed) == 1
                # return control_settings
                # map(x -> x.vals = abs.(round.(x.vals)), control_settings)
                # time_elapsed += @elapsed result = feasibility_pump_at_last(wm, control_settings)
                # control_settings = get_control_settings_from_result(result)
                # map(x -> x.vals = abs.(round.(x.vals)), control_settings)
                num_look_back += 2
            else
                num_look_back += 2
            end
        end
    end

    WM.Memento.info(LOGGER, "[HEUR] Could not find a heuristic solution.")
    return control_settings
    #return nothing
end


# function objective_min_losses(wm::WM.AbstractWaterModel, nws::Vector{Int})
#     objective = JuMP.AffExpr(0.0)

#     for nw in nws
#         for tank in values(WM.ref(wm, nw, :tank))
#             # head_min = WM.ref(wm, nw, :node, tank["node"], "head_min")
#             # head_max = WM.ref(wm, nw, :node, tank["node"], "head_max")
#             # head_mid = 0.5 * (head_min + head_max)
#             # objective += WM.var(wm, nw, :V, tank["index"])
#         end

#         # for pipe in values(WM.ref(wm, nw, :pipe))
#         #     dhp = WM.var(wm, nw, :dhp_pipe, pipe["index"])
#         #     dhn = WM.var(wm, nw, :dhn_pipe, pipe["index"])
#         #     objective += dhp + dhn
#         # end

#         for pump in values(WM.ref(wm, nw, :pump))
#             objective -= WM.var(wm, nw, :g_pump, pump["index"])
#             # objective += WM.var(wm, nw, :z_pump, pump["index"])
#         end

#         # for tank in values(WM.ref(wm, nw, :tank))
#         #     node = WM.ref(wm, nw, :node, tank["node"])
#         #     h = WM.var(wm, nw, :h, tank["node"])            
#         #     head_mid = node["head_min"] + 0.5 * (node["head_max"] - node["head_min"])
#         #     objective += (h - head_mid)^2
#         # end
#     end

#     JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective)
# end


function run_heuristic_mn(network_mn::Dict{String, <:Any}, mip_optimizer, nlp_optimizer, time_limit::Float64)
    network = WM.make_single_network(network_mn)
    WM.set_flow_partitions_si!(network_mn, 5.0, 1.0e-6)
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)

    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf)
    WM.JuMP.set_optimizer(wm.model, mip_optimizer)

    objective_mean_tank_levels(wm, sort(collect(WM.nw_ids(wm))))
    wm_cq = _instantiate_cq_model(network, nlp_optimizer)
    control_settings = run_sequential_heuristic(wm, wm_cq, result_micp, time_limit)

    if control_settings !== nothing
        results = simulate_control_settings_sequential(wm_cq, control_settings)
        is_feasible = all(x -> x.feasible, results)
        total_cost = sum(x.cost for x in results)
        WM.Memento.info(LOGGER, "[HEUR] Found heuristic solution, feasibility $(is_feasible).") 
        WM.Memento.info(LOGGER, "[HEUR] Found heuristic solution, cost $(total_cost).")
    end

    return control_settings
end


function run_heuristic(
    network_mn::Dict{String, <:Any}, network::Dict{String, <:Any},
    pc_path::String, mip_optimizer, nlp_optimizer, time_limit::Float64)
    WM.set_breakpoints!(network, 0.25, 1.0e-4)
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)

    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf)
    WM.JuMP.set_optimizer(wm.model, mip_optimizer)
    pairwise_cuts = load_pairwise_cuts(pc_path)
    add_pairwise_cuts(wm, pairwise_cuts)

    # objective_min_losses(wm, sort(collect(WM.nw_ids(wm)))[1:end-1])
    objective_mean_tank_levels(wm, sort(collect(WM.nw_ids(wm))))
    add_pump_volume_cuts!(wm)

    wm_cq = _instantiate_cq_model(network, nlp_optimizer)
    control_settings = run_sequential_heuristic(wm, wm_cq, result_micp, time_limit)
    
    if control_settings !== nothing
        results = simulate_control_settings_sequential(wm_cq, control_settings)
        is_feasible = all(x -> x.feasible, results)
        total_cost = sum(x.cost for x in results)
        WM.Memento.info(LOGGER, "[HEUR] Found heuristic solution, feasibility $(is_feasible).") 
        WM.Memento.info(LOGGER, "[HEUR] Found heuristic solution, cost $(total_cost).")
    end

    return control_settings
end


function run_heuristic_micp(network::Dict{String, Any}, mip_optimizer, nlp_optimizer, time_limit::Float64)
    network_mn = WM.make_multinetwork(network)
    wm_micp = construct_owf_model_micp(network_mn, nlp_optimizer, mip_optimizer)
    objective_min_losses(wm_micp, sort(collect(WM.nw_ids(wm_micp)))[1:end-1])
    return WM.optimize_model!(wm_micp; relax_integrality = false)
end


function run_heuristic(network::Dict{String, Any}, mip_optimizer, nlp_optimizer, time_limit::Float64)
    WM.set_breakpoints!(network, 0.25, 1.0e-4)
    network_mn = WM.make_multinetwork(network)
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)

    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf)
    WM.JuMP.set_optimizer(wm.model, mip_optimizer)

    # objective_min_losses(wm, sort(collect(WM.nw_ids(wm)))[1:end-1])
    objective_mean_tank_levels(wm, sort(collect(WM.nw_ids(wm))))
    add_pump_volume_cuts!(wm)
    # constraint_tank_final_range(wm)

    wm_cq = _instantiate_cq_model(network, nlp_optimizer)
    control_settings = run_sequential_heuristic(wm, wm_cq, result_micp, time_limit)

    if control_settings !== nothing
        results = simulate_control_settings_sequential(wm_cq, control_settings)
        is_feasible = all(x -> x.feasible, results)
        total_cost = sum(x.cost for x in results)
        WM.Memento.info(LOGGER, "[HEUR] Found heuristic solution, feasibility $(is_feasible).") 
        WM.Memento.info(LOGGER, "[HEUR] Found heuristic solution, cost $(total_cost).")
    end

    repair_schedule(control_settings, network, nlp_optimizer)
    return control_settings
end