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
    return ControlSetting(nw, vars, vals)
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

    one_vars = vars[findall(x -> round(x) == 1.0, vals)]
    zero_vars = vars[findall(x -> round(x) == 0.0, vals)]
    JuMP.@constraint(wm.model, sum(zero_vars) -
        sum(one_vars) >= 1.0 - length(one_vars))
end


function objective_mean_tank_levels(wm::WM.AbstractWaterModel)
    objective = JuMP.AffExpr(0.0)

    for nw in sort(collect(WM.nw_ids(wm)))
        for tank in values(WM.ref(wm, nw, :tank))
            node = WM.ref(wm, nw, :node, tank["node"])
            h = WM.var(wm, nw, :h, tank["node"])
            head_mid = 0.5 * (node["head_min"] + node["head_max"])
            objective += (h - head_mid)^2
        end
    end

    JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective)
end


function solve_pseudo_relaxation_nw(wm::WM.AbstractWaterModel, nws::Vector{Int})
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    nw_last_binary = sort(findall(x -> x in nws, nw_ids))[end]
    nw_relaxed = nw_ids[nw_last_binary+1:end]

    vars_to_relax = vcat(WM.get_all_binary_vars_at_nw!.(Ref(wm), nw_relaxed)...)
    map(x -> !JuMP.is_fixed(x) && JuMP.unset_binary(x), vars_to_relax)
    map(x -> !JuMP.is_fixed(x) && JuMP.set_lower_bound(x, 0.0), vars_to_relax)
    map(x -> !JuMP.is_fixed(x) && JuMP.set_upper_bound(x, 1.0), vars_to_relax)

    result = WM.optimize_model!(wm; relax_integrality = false)
    map(x -> !JuMP.is_binary(x) && JuMP.set_binary(x), vars_to_relax)

    return result
end


function run_sequential_heuristic(
    wm::WM.AbstractWaterModel, wm_cq::CQWaterModel,
    result_micp::Dict{String, <:Any}, time_limit::Float64)
    control_settings = get_control_settings_from_result(result_micp)
    map(x -> x.vals = round.(x.vals), control_settings)
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    time_elapsed = 0.0

    while time_elapsed < time_limit
        # Simulate the current control schedule forward in time.
        time_elapsed += @elapsed simulation_results =
            simulate_control_settings_sequential(wm_cq, control_settings)

        # If the control schedule is feasible...
        if all(x -> x.feasible, simulation_results)
            # Return the control settings schedule.
            return control_settings
        else # If the control schedule is not feasible...
            # Get the set of times at and before the infeasibility.
            nws_infeasible = nw_ids[1:length(simulation_results)]
            WM.Memento.info(LOGGER, "[HEUR] Feasible time indices: $(nws_infeasible).")

            # Add a feasibility cut to the WaterModels model up to this time.
            add_feasibility_cut_from_control_settings!(
                wm, control_settings[nws_infeasible])

            # Pick a random time index at which the controls will be updated.
            num_choose = min(length(nws_infeasible), 1)
            nws_random = sort(Random.shuffle(nws_infeasible)[1:num_choose])

            # Unfix the controls at this time step and beyond.
            nws_unfixed = nw_ids[nws_random[end]:end]
            unfix_control_setting!.(Ref(wm), control_settings[nws_unfixed])
 
            # Fix the controls prior to this time step.
            nws_fixed = nw_ids[1:nws_random[end]-1]
            fix_control_setting!.(Ref(wm), control_settings[nws_fixed],
                simulation_results[nws_fixed])

            # Solve the corresponding pseudo-relaxation at the random time.
            time_elapsed += @elapsed result =
                solve_pseudo_relaxation_nw(wm, nws_random)

            if result["primal_status"] === WM._MOI.FEASIBLE_POINT
                control_settings = get_control_settings_from_result(result)
                map(x -> x.vals = round.(x.vals), control_settings)
            else
                add_feasibility_cut_from_control_settings!(
                    wm, control_settings[nws_infeasible])
            end
        end
    end

    WM.Memento.info(LOGGER, "[HEUR] Could not find a heuristic solution.")
    return nothing
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

    objective_mean_tank_levels(wm)
    add_pump_volume_cuts!(wm)

    wm_cq = _instantiate_cq_model(network, nlp_optimizer)
    return run_sequential_heuristic(wm, wm_cq, result_micp, time_limit)
end


function run_heuristic(network::Dict{String, Any}, mip_optimizer, nlp_optimizer, time_limit::Float64)
    set_breakpoints_accuracy!(network, 0.25, 1.0e-4)
    network_mn = WM.make_multinetwork(network)
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)

    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf)
    WM.JuMP.set_optimizer(wm.model, mip_optimizer)

    objective_mean_tank_levels(wm)
    add_pump_volume_cuts!(wm)

    wm_cq = _instantiate_cq_model(network, nlp_optimizer)
    return run_sequential_heuristic(wm, wm_cq, result_micp, time_limit)

    # # Specify model options and construct the multinetwork OWF model.
    # wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf)

    # # Set the optimizer and other important solver parameters.
    # WM.JuMP.set_optimizer(wm.model, mip_optimizer)

    # binary_vars = filter(v -> JuMP.is_binary(v), JuMP.all_variables(wm.model))
    # vars_to_relax_not_y = filter(v -> !occursin("_y", JuMP.name(v)), binary_vars)
    # vars_to_relax_not_xy = filter(v -> !occursin("_x", JuMP.name(v)), vars_to_relax_not_y)
    # map(x -> !JuMP.is_fixed(x) && JuMP.unset_binary(x), vars_to_relax_not_xy)

    # return WM.optimize_model!(wm)


    # set_breakpoints_num_mn!(network_mn, 10)

    # networks_mn = [deepcopy(network_mn) for i in 1:Threads.nthreads()]
    # networks = [deepcopy(network) for i in 1:Threads.nthreads()]
    # heuristic_result = Vector{Any}([nothing for i in 1:Threads.nthreads()])
    # # wms = [_instantiate_cq_model(networks[i], nlp_optimizer) for i in 1:Threads.nthreads()]
    # # settings = [[ControlSetting(n, [], [])] for n in 1:Threads.nthreads()]
    # # results = [[SimulationResult(false, Dict{Int, Float64}(), 0.0)] for n in 1:Threads.nthreads()]
    # feasible, costs = [false for n in 1:Threads.nthreads()], zeros(Float64, Threads.nthreads())
    
    # Threads.@threads for k in 1:Threads.nthreads()
    #     heuristic_result[k] = run_feasibility_pump(networks_mn[k], mip_optimizer)
    #     is_feasible = heuristic_result[k] === nothing ? false : 
    #         heuristic_result[k]["primal_status"] === WM._MOI.FEASIBLE_POINT

    #     println(is_feasible)

    #     if is_feasible
    #         feasible[k], costs[k] = WaterModelsAnnex.simulate_result_mn(
    #             networks[k], heuristic_result[k], nlp_optimizer)
    #     else
    #         feasible[k], costs[k] = false, Inf
    #     end

    #     # settings[k] = is_feasible ? get_control_settings_from_result(heuristic_result[k]) : settings[k]
    #     # results[k] = is_feasible ? simulate_control_setting.(Ref(wms[k]), settings[k]) : results[k]
    #     # feasible[k] = is_feasible ? all(x.feasible for x in results[k]) : feasible[k]
    #     # costs[k] = is_feasible ? sum(x.cost for x in results[k]) : Inf
    #     # costs[k] = is_feasible ? heuristic_result[k]["objective"] : Inf
    # end

    # println("Heuristic solution costs: ", minimum(costs), " ", costs)
    # heuristic_result[findmin(costs)[2]]["objective"] = minimum(costs)
    # return heuristic_result[findmin(costs)[2]]
end