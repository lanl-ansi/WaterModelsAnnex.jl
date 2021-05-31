function solve_obbt(network_path::String, modification_path::String, obbt_optimizer; time_limit::Float64 = 3600.0, use_obbt::Bool = true, kwargs...)
    # Read in the original network data.
    network = WM.parse_file(network_path; skip_correct = true)
    modifications = WM.parse_file(modification_path; skip_correct = true)
    WM._IM.update_data!(network, modifications)
    WM.correct_network_data!(network)
    num_obbt_rounds = use_obbt ? 100 : 1

    # Tighten bounds of variables in the network.
    WM.solve_obbt_owf!(
        network, obbt_optimizer; model_type = WM.PWLRDWaterModel,
        solve_relaxed = false, time_limit = time_limit, max_iter = num_obbt_rounds,
        ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10))

    # Get tightened network data.
    return network
end


function solve_owf(network_path::String, modification_path::String, obbt_optimizer, owf_optimizer, nlp_optimizer; kwargs...)
    # Print the amount of time spent to perform the above.
    WM.Memento.info(LOGGER, "Starting optimization-based bound tightening.")

    # Tighten bounds in the network using optimization-based bound tightening.
    obbt_time = @elapsed network = solve_obbt(network_path,
        modification_path, obbt_optimizer; kwargs...)

    # Print the amount of time spent to perform the above.
    WM.Memento.info(LOGGER, "Optimization-based bound tightening completed in $(obbt_time) seconds.")

    # Solve the OWF problem using the network from above.
    return solve_owf(network_path, network, obbt_optimizer,
        owf_optimizer, nlp_optimizer; kwargs...)
end


function compute_pairwise_cuts(network::Dict{String, Any}, optimizer)
    # Relax the network.
    WM._relax_network!(network)

    # Specify extensions to be used in the WaterModels formulation.
    ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10)

    # Construct independent thread-local WaterModels objects.
    wms = [WM.instantiate_model(network, WM.PWLRDWaterModel,
        WM.build_owf; ext = ext) for i in 1:Threads.nthreads()]

    # Set the optimizer for the WaterModels objects.
    map(x -> JuMP.set_optimizer(x.model, optimizer), wms)

    # Get problem sets for generating pairwise cuts.
    problem_sets = WM._get_pairwise_problem_sets(wms[1])

    # Generate data structures to store thread-local cut results.
    cuts_array = Array{Array{WM._PairwiseCut, 1}, 1}([])

    for i in 1:Threads.nthreads()
        # Initialize the per-thread pairwise cut array.
        push!(cuts_array, Array{WM._PairwiseCut, 1}([]))
    end

    Threads.@threads for i in 1:length(problem_sets)
        # Compute pairwise cuts across all problem sets.
        cuts_local = WM._compute_pairwise_cuts!(wms[Threads.threadid()], [problem_sets[i]])
        append!(cuts_array[Threads.threadid()], cuts_local)
    end

    # Concatenate all cutting plane results.
    cuts = vcat(cuts_array...)

    # Unrelax the network.
    WM._fix_demands!(network)
    WM._fix_tanks!(network)
    WM._fix_reservoirs!(network)

    # Return the data structure comprising cuts.
    return cuts
end


function construct_owf_model_relaxed(network_mn::Dict{String, Any}, optimizer; kwargs...)
    wm = WM.instantiate_model(network_mn, WM.CRDWaterModel, WM.build_mn_owf)
    WM.JuMP.set_optimizer(wm.model, optimizer)
    return wm # Return the relaxation-based WaterModels object.
end


function construct_owf_model(network_mn::Dict{String, Any}, owf_optimizer; use_new::Bool = true, kwargs...)
    # Specify model options and construct the multinetwork OWF model.
    model_type = use_new ? WM.PWLRDWaterModel : WM.LRDWaterModel
    wm = WM.instantiate_model(network_mn, model_type, WM.build_mn_owf)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, owf_optimizer)
    WM._MOI.set(wm.model, WM._MOI.NumberOfThreads(), 1)

    # Return the WaterModels object.
    return wm
end


function add_pairwise_cuts(wm::WM.AbstractWaterModel, cuts::Array{WM._PairwiseCut, 1})
    # Add the pairwise cuts obtained from the relaxed problem to the OWF problem.
    for nw_id in sort(collect(WM.nw_ids(wm)))[1:end-1]
        # Use the same cuts for all subnetworks of the multinetwork.
        map(x -> x.variable_index_1.network_index = nw_id, cuts)
        map(x -> x.variable_index_2.network_index = nw_id, cuts)

        # Add the collection of pairwise cuts for the subnetwork.
        WM._add_pairwise_cuts!(wm, cuts)
    end
end


function degree_fr(data::Dict{String, <:Any}, component::Dict{String, <:Any})
    node_fr = component["node_fr"]

    pipes_fr = filter(x -> x.second["node_fr"] == node_fr, data["pipe"])
    pumps_fr = filter(x -> x.second["node_fr"] == node_fr, data["pump"])
    regulators_fr = filter(x -> x.second["node_fr"] == node_fr, data["regulator"])
    short_pipes_fr = filter(x -> x.second["node_fr"] == node_fr, data["short_pipe"])
    valves_fr = filter(x -> x.second["node_fr"] == node_fr, data["valve"])

    return length(pipes_fr) + length(pumps_fr) #+ length(regulators_fr) +
        #length(short_pipes_fr) + length(valves_fr)
end


function degree_to(data::Dict{String, <:Any}, component::Dict{String, <:Any})
    node_to = component["node_to"]

    pipes_to = filter(x -> x.second["node_to"] == node_to, data["pipe"])
    pumps_to = filter(x -> x.second["node_to"] == node_to, data["pump"])
    regulators_to = filter(x -> x.second["node_to"] == node_to, data["regulator"])
    short_pipes_to = filter(x -> x.second["node_to"] == node_to, data["short_pipe"])
    valves_to = filter(x -> x.second["node_to"] == node_to, data["valve"])

    return length(pipes_to) + length(pumps_to) #+ length(regulators_to) +
        #length(short_pipes_to) + length(valves_to)
end


function update_multinetwork_breakpoints!(network_mn::Dict{String, <:Any}, result_mn::Dict{String, <:Any})
    for nw_id in sort([parse(Int, x) for x in collect(keys(network_mn["nw"]))])[1:end-1]
        network = network_mn["nw"][string(nw_id)]
        result = result_mn["solution"]["nw"][string(nw_id)]

        for (i, pipe) in network["pipe"]
            flow_result = result["pipe"][i]["q"]
            flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]

            if degree_fr(network, pipe) + degree_to(network, pipe) <= 2
                pipe["flow_lower_breakpoints"] = [flow_min, flow_result, flow_max]
                pipe["flow_upper_breakpoints"] = [flow_min, flow_result, flow_max]
            else
                pipe["flow_lower_breakpoints"] = [flow_min, flow_max]
                pipe["flow_upper_breakpoints"] = [flow_min, flow_max]
            end
         end

         for (i, pump) in network["pump"]
            flow_result = max(0.0, result["pump"][i]["q"])
            flow_min, flow_max = pump["flow_min"], pump["flow_max"]

            if flow_result != 0.0
                pump["flow_lower_breakpoints"] = [flow_min, flow_result, flow_max]
                pump["flow_upper_breakpoints"] = [flow_min, flow_result, flow_max]
            else
                flow_mid = flow_min + 0.5 * (flow_max - flow_min)
                pump["flow_lower_breakpoints"] = [flow_min, flow_mid, flow_max]
                pump["flow_upper_breakpoints"] = [flow_min, flow_mid, flow_max]
            end
         end
    end
end


function update_multinetwork_heuristic_breakpoints!(network_mn::Dict{String, <:Any}, result_mn::Dict{String, <:Any})
    for nw_id in sort([parse(Int, x) for x in collect(keys(network_mn["nw"]))])[1:end-1]
        network = network_mn["nw"][string(nw_id)]
        result = result_mn["solution"]["nw"][string(nw_id)]

        for (i, pipe) in network["pipe"]
            flow_result = result["pipe"][i]["q"]
            flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]

            # pipe["flow_lower_breakpoints"] = range(flow_min, stop = flow_max, length = 5)
            # pipe["flow_upper_breakpoints"] = range(flow_min, stop = flow_max, length = 5)

            if degree_fr(network, pipe) + degree_to(network, pipe) <= 2
                flow_min_mid = flow_min + 0.5 * (flow_result - flow_min)
                flow_max_mid = flow_result + 0.5 * (flow_max - flow_result)
                pipe["flow_lower_breakpoints"] = [flow_min, flow_min_mid, flow_result, flow_max_mid, flow_max]
                pipe["flow_upper_breakpoints"] = [flow_min, flow_min_mid, flow_result, flow_max_mid, flow_max]
            else
                flow_mid = flow_min + 0.5 * (flow_max - flow_min)
                flow_min_mid = flow_min + 0.5 * (flow_mid - flow_min)
                flow_max_mid = flow_mid + 0.5 * (flow_max - flow_mid)
                pipe["flow_lower_breakpoints"] = [flow_min, flow_mid, flow_max]
                pipe["flow_upper_breakpoints"] = [flow_min, flow_max]
            end
         end

         for (i, pump) in network["pump"]
            flow_result = max(0.0, result["pump"][i]["q"])
            flow_min, flow_max = pump["flow_min_forward"], pump["flow_max"]

            # pump["flow_lower_breakpoints"] = range(flow_min, stop = flow_max, length = 5)
            # pump["flow_upper_breakpoints"] = range(flow_min, stop = flow_max, length = 5)

            if flow_result != 0.0
                flow_min_mid = flow_min + 0.5 * (flow_result - flow_min)
                flow_max_mid = flow_result + 0.5 * (flow_max - flow_result)
                pump["flow_lower_breakpoints"] = [flow_min, flow_min_mid, flow_result, flow_max_mid, flow_max]
                pump["flow_upper_breakpoints"] = [flow_min, flow_min_mid, flow_result, flow_max_mid, flow_max]
            else
                flow_mid = flow_min + 0.5 * (flow_max - flow_min)
                flow_min_mid = flow_min + 0.5 * (flow_mid - flow_min)
                flow_max_mid = flow_mid + 0.5 * (flow_max - flow_mid)
                pump["flow_lower_breakpoints"] = [flow_min, flow_max]
                pump["flow_upper_breakpoints"] = [flow_min, flow_mid, flow_max]
            end
         end
    end
end


function set_branching_priorities!(wm::WM.AbstractWaterModel)
    priority = 1

    for comp_type in [:pump, :regulator, :valve]
        for nw in sort(collect(WM.nw_ids(wm)))[1:end-1]
            for i in WM.ids(wm, nw, comp_type)
                var_symbol = Symbol("z_" * string(comp_type))
                var = WM.var(wm, nw, var_symbol, i)

                WM._MOI.set(JuMP.backend(wm.model).optimizer,
                    Gurobi.VariableAttribute("BranchPriority"),
                    JuMP.index(var), priority)
            end

            priority += 1
        end
    end

    for comp_type in [:pipe, :pump, :regulator, :short_pipe, :valve]
        for nw in sort(collect(WM.nw_ids(wm)))[1:end-1]
            for i in WM.ids(wm, nw, comp_type)
                var_symbol = Symbol("y_" * string(comp_type))
                var = WM.var(wm, nw, var_symbol, i)
                
                WM._MOI.set(JuMP.backend(wm.model).optimizer,
                    Gurobi.VariableAttribute("BranchPriority"),
                    JuMP.index(var), priority)
            end

            priority += 1
        end
    end
end


function solve_owf(network_path::String, network, obbt_optimizer, owf_optimizer, nlp_optimizer; use_new::Bool = true, kwargs...)
    network_mn = WM.make_multinetwork(network)
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)
    update_multinetwork_breakpoints!(network_mn, result_micp)

    # Construct the OWF model that will serve as the master problem.
    wm_relaxed = construct_owf_model(network_mn, obbt_optimizer; use_new = use_new, kwargs...)
    result_relaxed = WM.optimize_model!(wm_relaxed; relax_integrality = true)

    _update_tank_time_series_heur!(network, result_micp)
    wm_master = construct_owf_model(network_mn, owf_optimizer; use_new = use_new, kwargs...)
    heuristic_setting = calc_heuristic(wm_master, network, obbt_optimizer, nlp_optimizer)
    heur_cost = set_warm_start_from_setting!(wm_master, network, heuristic_setting, nlp_optimizer)
    WM.Memento.info(LOGGER, "Heuristic found solution with objective $(round(heur_cost; digits = 2)).")

    if use_new
        # Add binary-binary and binary-continuous pairwise cuts.
        WM.Memento.info(LOGGER, "Beginning cut preprocessing routine.")
        cut_time = @elapsed pairwise_cuts = compute_pairwise_cuts(network, obbt_optimizer)
        cut_time += @elapsed add_pairwise_cuts(wm_master, pairwise_cuts)
        cut_time += @elapsed add_pump_volume_cuts!(wm_master)
        WM.Memento.info(LOGGER, "Cut preprocessing completed in $(cut_time) seconds.")
    end

    # TODO: Remove this once Gurobi.jl interface is fixed.
    wm_master.model.moi_backend.optimizer.model.has_generic_callback = false

    # Add the lazy cut callback.
    lazy_cut_stats = add_owf_lazy_cut_callback!(wm_master, network, heuristic_setting[1], nlp_optimizer)

    if use_new
        # Add the user cut callback.
        WM.Memento.info(LOGGER, "Building and applying user cut callback routine.")
        user_cut_stats = add_owf_user_cut_callback!(wm_master)
    end

    # Add the heuristic callback.
    # add_owf_heuristic_callback!(wm_master, network, nlp_optimizer, obbt_optimizer)

    # Optimize the master WaterModels model.
    WM.Memento.info(LOGGER, "Solving the master OWF problem.")
    WM._relax_all_direction_variables!(wm_master)

    result = WM.optimize_model!(wm_master; relax_integrality = false)
    WM.Memento.info(LOGGER, "Solved for $(result["solve_time"]) seconds.")

    set_branching_priorities!(wm_master)
    result = WM.optimize_model!(wm_master; relax_integrality = false)
    WM.Memento.info(LOGGER, "Solved for $(result["solve_time"]) seconds.")

    # Print out relevant statistics of the solution process.
    percent_lazy = round(lazy_cut_stats.time_elapsed / result["solve_time"] * 100.0; digits = 2)
    WM.Memento.info(LOGGER, "Lazy cuts accounted for $(percent_lazy)% of solve time.")

    if use_new
        percent_user = round(user_cut_stats.time_elapsed / result["solve_time"] * 100.0; digits = 2)
        WM.Memento.info(LOGGER, "User cuts accounted for $(percent_user)% of solve time.")
    end

    return result
end


function filter_feasible_control_settings!(control_settings::Array{ControlSetting, 1}, control_results::Array{SimulationResult, 1})
    feasible_indices = findall(x -> x.feasible, control_results)
    control_settings = control_settings[feasible_indices]
    control_results = control_results[feasible_indices]
end


function control_settings_from_solution(control_sol::Dict{Int, <:Any}, settings::Array{ControlSetting, 1})
    settings = [deepcopy(settings[1]) for i in 1:length(keys(control_sol))]

    for (k, nw) in enumerate(sort(collect(keys(control_sol))))
        settings[k].network_id = nw
        settings[k].vals = abs.(control_sol[nw])
    end

    return settings
end


function calc_heuristic(wm::WM.AbstractWaterModel, network::Dict{String, <:Any}, mip_optimizer, nlp_optimizer)
    # Build all possible control settings, simulate them, and filter the results.
    WM.Memento.info(LOGGER, "Precomputing all possible heuristic control settings.")
    settings = create_all_control_settings(wm)

    network_ids = sort(unique([x.network_id for x in settings]))
    return [deepcopy(settings[1]) for i in 1:length(network_ids)]

    WM.Memento.info(LOGGER, "Simulating all heuristic control settings.")
    results = simulate_control_settings(network, settings, nlp_optimizer)

    WM.Memento.info(LOGGER, "Filtering feasible heuristic control settings.")
    filter_feasible_control_settings!(settings, results)

    WM.Memento.info(LOGGER, "Solving heuristic linear programming relaxation.")
    weights = solve_heuristic_linear_program(wm, settings, results, nlp_optimizer)

    WM.Memento.info(LOGGER, "Solving heuristic master mixed-integer program.")
    control_sol = solve_heuristic_master_program(wm, network, settings, weights, mip_optimizer, nlp_optimizer)

    if control_sol !== nothing
        return control_settings_from_solution(control_sol, settings)
    else
        network_ids = sort(unique([x.network_id for x in settings]))
        return [deepcopy(settings[1]) for i in 1:length(network_ids)]
    end
end
