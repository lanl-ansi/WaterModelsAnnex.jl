function compute_pairwise_cuts(network::Dict{String, Any}, optimizer)
    # Relax the network.
    WM._relax_network!(network)

    # Construct independent thread-local WaterModels objects.
    wms = [WM.instantiate_model(network, WM.PWLRDWaterModel,
        WM.build_owf) for i in 1:Threads.nthreads()]

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


# function construct_owf_model_micp(network_mn::Dict{String, Any}, nlp_optimizer, mip_optimizer; kwargs...)
#     wm = WM.instantiate_model(network_mn, WM.CRDWaterModel, WM.build_mn_owf)
#     f = Juniper.register(WM.head_loss_args(wm)..., autodiff=false)
#     optimizer = JuMP.optimizer_with_attributes(Juniper.Optimizer,
#             "nl_solver" => nlp_optimizer, "registered_functions" => [f],
#             "allow_almost_solved_integral" => false)
#     WM.JuMP.set_optimizer(wm.model, optimizer)
#     return wm # Return the relaxation-based WaterModels object.
# end


function construct_owf_model(network_mn::Dict{String, Any}, owf_optimizer; use_pwlrd::Bool = true, kwargs...)
    # Specify model options and construct the multinetwork OWF model.
    model_type = use_pwlrd ? WM.PWLRDWaterModel : LRDXWaterModel
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


function set_lazy_attributes!(wm::WM.AbstractWaterModel)
    for nw in sort(collect(WM.nw_ids(wm)))[1:end-1]
        for constraints in values(WM.con(wm, nw, :pipe_head_loss))
            WM._MOI.set.(Ref(JuMP.backend(wm.model).optimizer),
                    Ref(Gurobi.ConstraintAttribute("Lazy")),
                    JuMP.index.(constraints), 2)
        end

        for constraints in values(WM.con(wm, nw, :on_off_pump_head_gain))
            WM._MOI.set.(Ref(JuMP.backend(wm.model).optimizer),
                    Ref(Gurobi.ConstraintAttribute("Lazy")),
                    JuMP.index.(constraints), 2)
        end
    end
end


function set_branching_priorities!(wm::WM.AbstractWaterModel)
    priority = 1

    for comp_type in [:valve, :short_pipe, :regulator, :pump, :pipe]
        for nw in reverse(sort(collect(WM.nw_ids(wm)))[1:end-1])
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

    for comp_type in [:valve, :regulator, :pump]
        for nw in reverse(sort(collect(WM.nw_ids(wm)))[1:end-1])
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
end


function solve_owf(network::Dict, mip_optimizer, nlp_optimizer, breakpoint_function!::Function)
    # Parse the network data.
    network_mn = WM.make_multinetwork(network)

    # Solve a continuously-relaxed version of the problem.
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)

    # Set the breakpoints to be used for nonlinear functions.
    breakpoint_function!(network_mn, result_micp)
    wm_master = construct_owf_model(network_mn, mip_optimizer)

    # Solve the model and return the result.
    return WM.optimize_model!(wm_master; relax_integrality = false)
end


function solve_owf_with_cuts(network::Dict, pc_path::String, mip_optimizer, nlp_optimizer, breakpoint_function!::Function, cut_types::Int)
    # Parse the network data.
    network_mn = WM.make_multinetwork(network)

    # Solve a continuously-relaxed version of the problem.
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)

    # Set the breakpoints to be used for nonlinear functions.
    breakpoint_function!(network_mn, result_micp)
    wm_master = construct_owf_model(network_mn, mip_optimizer)

    if cut_types == 1
        pairwise_cuts = load_pairwise_cuts(pc_path)
        add_pairwise_cuts(wm_master, pairwise_cuts)
    elseif cut_types == 2
        pairwise_cuts = load_pairwise_cuts(pc_path)
        add_pairwise_cuts(wm_master, pairwise_cuts)
        add_pump_volume_cuts!(wm_master)
    end

    # Solve the model and return the result.
    return WM.optimize_model!(wm_master; relax_integrality = false)
end


function solve_owf_formulation(network::Dict, pc_path::String, mip_optimizer, nlp_optimizer, breakpoint_function!::Function; use_pwlrd::Bool = true)
    # Parse the network data.
    network_mn = WM.make_multinetwork(network)

    # Solve a continuously-relaxed version of the problem.
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)

    # Set the breakpoints to be used for nonlinear functions.
    breakpoint_function!(network_mn, result_micp)
    wm_master = construct_owf_model(network_mn, mip_optimizer; use_pwlrd = use_pwlrd)

    pairwise_cuts = load_pairwise_cuts(pc_path)
    add_pairwise_cuts(wm_master, pairwise_cuts)
    add_pump_volume_cuts!(wm_master)

    # Solve the model and return the result.
    return WM.optimize_model!(wm_master; relax_integrality = false)
end


function solve_owf_upper_bounds(network_mn::Dict, pc_path::String, mip_optimizer, nlp_optimizer)
    # Solve a continuously-relaxed version of the problem.
    wm_micp = WM.instantiate_model(network_mn, WM.CRDWaterModel, WM.build_mn_owf)
    WM.JuMP.set_optimizer(wm_micp.model, nlp_optimizer)
    WM.optimize_model!(wm_micp; relax_integrality = true)

    # Set the optimizer and other important solver parameters.
    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf)
    WM.JuMP.set_optimizer(wm.model, mip_optimizer)
    WM._MOI.set(wm.model, WM._MOI.NumberOfThreads(), 1)

    # Add extra cuts to the formulation.
    pairwise_cuts = load_pairwise_cuts(pc_path)
    add_pairwise_cuts(wm, pairwise_cuts)
    add_pump_volume_cuts!(wm)

    # TODO: Remove this once Gurobi.jl interface is fixed.
    wm.model.moi_backend.optimizer.model.has_generic_callback = false

    # Setup a single-step version of the network.
    network = WM.make_single_network(network_mn)

    # Add the lazy cut callback.
    lazy_cut_stats = add_owf_lazy_cut_callback!(wm, network, nlp_optimizer)

    # Solve the model and return the result.
    result = WM.optimize_model!(wm; relax_integrality = false)
    result["true_upper_bound"] = lazy_cut_stats.best_cost
    result["true_gap"] = (result["true_upper_bound"] -
        result["objective_lb"]) / result["true_upper_bound"]
    return result
end


function solve_owf_upper_bounds(network_mn::Dict, network::Dict, build_method::Function, mip_optimizer, nlp_optimizer)
    # Solve a continuously-relaxed version of the problem.
    wm_micp = WM.instantiate_model(network_mn, WM.CRDWaterModel, build_method)
    WM.JuMP.set_optimizer(wm_micp.model, nlp_optimizer)
    WM.optimize_model!(wm_micp; relax_integrality = true)
 
    # Set the optimizer and other important solver parameters.
    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, build_method)
    WM.JuMP.set_optimizer(wm.model, mip_optimizer)
    WM._MOI.set(wm.model, WM._MOI.NumberOfThreads(), 1)

    # TODO: Remove this once Gurobi.jl interface is fixed.
    wm.model.moi_backend.optimizer.model.has_generic_callback = false

    # Add the lazy cut callback.
    lazy_cut_stats = add_owf_lazy_cut_callback!(wm, network, nlp_optimizer)    

    # Solve the model and store the result.
    result = WM.optimize_model!(wm)

    # Add true upper bound data to the result dictionary.
    result["true_upper_bound"] = lazy_cut_stats.best_cost
    result["true_gap"] = (result["true_upper_bound"] -
        result["objective_lb"]) / result["true_upper_bound"]

    # Return the result dictionary.
    return result
end


function solve_owf_upper_bounds_after_repair(network::Dict, mip_optimizer, nlp_optimizer)
    # Parse the network data.
    network_mn = WM.make_multinetwork(network)

    # Solve a continuously-relaxed version of the problem.
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)
    control_settings = get_control_settings_from_result(result_micp)
    repair_schedule(control_settings, network, nlp_optimizer)

    wm_master = construct_owf_model(network_mn, mip_optimizer; use_pwlrd = true)
    set_warm_start_from_setting!(wm_master, network, control_settings, nlp_optimizer)

    wm_master.model.moi_backend.optimizer.model.has_generic_callback = false
    lazy_cut_stats = add_owf_lazy_cut_callback!(wm_master, network, control_settings[1], nlp_optimizer)
    result = WM.optimize_model!(wm_master; relax_integrality = false)
    
    
    # direction_settings = get_direction_settings_from_result(result_micp)
    # # Set the breakpoints to be used for nonlinear functions.
    # if formulation_type == 1
    #     set_breakpoints_piecewise_degree!(network_mn, result_micp)
    #     wm_master = construct_owf_model(network_mn, mip_optimizer; use_pwlrd = false)
    #     # pairwise_cuts = load_pairwise_cuts(pc_path)
    #     # add_pairwise_cuts(wm_master, pairwise_cuts)
    #     add_pump_volume_cuts!(wm_master)
    # elseif formulation_type == 2
    #     set_breakpoints_oa!(network_mn, result_micp)
    #     wm_master = construct_owf_model(network_mn, mip_optimizer; use_pwlrd = true)
    #     # WM._relax_all_direction_variables!(wm_master)
    # end

    # # TODO: Remove this once Gurobi.jl interface is fixed.
    # wm_master.model.moi_backend.optimizer.model.has_generic_callback = false

    # # Add the lazy cut callback.
    # lazy_cut_stats = add_owf_lazy_cut_callback!(wm_master, network,
    #     control_settings[1], nlp_optimizer)    

    # # Solve the model and return the result.
    # result = WM.optimize_model!(wm_master; relax_integrality = false)
    # result["true_upper_bound"] = lazy_cut_stats.best_cost
    # result["true_gap"] = (result["true_upper_bound"] - result["objective_lb"]) / result["true_upper_bound"]
    # return result
end


function solve_owf_upper_bounds_bp(network::Dict, pc_path::String, mip_optimizer, nlp_optimizer, error_tolerance::Float64)
    # Parse the network data.
    WM.set_breakpoints!(network, error_tolerance, 1.0e-4)
    network_mn = WM.make_multinetwork(network)

    # Solve a continuously-relaxed version of the problem.
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)
    control_settings = get_control_settings_from_result(result_micp)

    # Set the breakpoints to be used for nonlinear functions.
    wm_master = construct_owf_model(network_mn, mip_optimizer; use_pwlrd = false)
    pairwise_cuts = load_pairwise_cuts(pc_path)
    add_pairwise_cuts(wm_master, pairwise_cuts)
    add_pump_volume_cuts!(wm_master)

    # TODO: Remove this once Gurobi.jl interface is fixed.
    wm_master.model.moi_backend.optimizer.model.has_generic_callback = false

    # Add the lazy cut callback.
    lazy_cut_stats = add_owf_lazy_cut_callback!(
        wm_master, network, control_settings[1], nlp_optimizer)    

    # Solve the model and return the result.
    result = WM.optimize_model!(wm_master; relax_integrality = false)
    result["true_upper_bound"] = lazy_cut_stats.best_cost
    result["true_gap"] = (result["true_upper_bound"] -
        result["objective_lb"]) / result["true_upper_bound"]

    return result
end


function solve_owf_upper_bounds_no_cb(network::Dict, pc_path::String, mip_optimizer, nlp_optimizer, formulation_type::Int)
    # Parse the network data.
    network_mn = WM.make_multinetwork(network)

    # Solve a continuously-relaxed version of the problem.
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)
    control_settings = get_control_settings_from_result(result_micp)

    # Set the breakpoints to be used for nonlinear functions.
    if formulation_type == 1
        set_breakpoints_piecewise_degree!(network_mn, result_micp)
        wm_master = construct_owf_model(network_mn, mip_optimizer; use_pwlrd = false)
        pairwise_cuts = load_pairwise_cuts(pc_path)
        add_pairwise_cuts(wm_master, pairwise_cuts)
        add_pump_volume_cuts!(wm_master)
    elseif formulation_type == 2
        set_breakpoints_oa!(network_mn, result_micp)
        wm_master = construct_owf_model(network_mn, mip_optimizer; use_pwlrd = true)
        WM._relax_all_direction_variables!(wm_master)
    end

    # Solve the model and return the result.
    return WM.optimize_model!(wm_master; relax_integrality = false)
end


function solve_pairwise_cuts(network::Dict, optimizer)
    WM.Memento.info(LOGGER, "Beginning cut preprocessing routine.")
    cut_time = @elapsed pairwise_cuts = compute_pairwise_cuts(network, optimizer)
    WM.Memento.info(LOGGER, "Pairwise cut preprocessing completed in $(cut_time) seconds.")
    return pairwise_cuts
end


function load_pairwise_cuts(path::String)
    cuts_array = Vector{WM._PairwiseCut}([])
    
    for entry in WM.JSON.parsefile(path)
        vid_1_network_index = Int(entry["variable_index_1"]["network_index"])
        vid_1_component_type = Symbol(entry["variable_index_1"]["component_type"])
        vid_1_variable_symbol = Symbol(entry["variable_index_1"]["variable_symbol"])
        vid_1_component_index = Int(entry["variable_index_1"]["component_index"])
        vid_1 = WM._VariableIndex(vid_1_network_index, vid_1_component_type,
            vid_1_variable_symbol, vid_1_component_index)

        vid_2_network_index = Int(entry["variable_index_2"]["network_index"])
        vid_2_component_type = Symbol(entry["variable_index_2"]["component_type"])
        vid_2_variable_symbol = Symbol(entry["variable_index_2"]["variable_symbol"])
        vid_2_component_index = Int(entry["variable_index_2"]["component_index"])
        vid_2 = WM._VariableIndex(vid_2_network_index, vid_2_component_type,
            vid_2_variable_symbol, vid_2_component_index)

        coefficient_1 = entry["coefficient_1"]
        coefficient_2 = entry["coefficient_2"]
        constant = entry["constant"]
        
        push!(cuts_array, WM._PairwiseCut(coefficient_1, vid_1,
            coefficient_2, vid_2, constant))
    end

    return cuts_array
end


function simulate_result_mn(network::Dict{String, <:Any}, result::Dict{String, <:Any}, nlp_optimizer)
    wm_cq = _instantiate_cq_model(deepcopy(network), nlp_optimizer)
    control_settings = get_control_settings_from_result(deepcopy(result))
    simulate_control_settings_sequential(wm_cq, control_settings)
end


function solve_owf_with_heuristic(network::Dict, heuristic_path::String, mip_optimizer, nlp_optimizer)
    # Solve a relaxed version of the problem and set breakpoints.
    network_mn = WM.make_multinetwork(network)
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)
    set_breakpoints_piecewise_degree!(network_mn, result_micp)

    heuristic_result = WM.JSON.parsefile(heuristic_path)
    wm_cq = _instantiate_cq_model(network, nlp_optimizer)
    control_settings = get_control_settings_from_result(heuristic_result)

    # Solve the model and return the result.
    # wm_master = construct_owf_model(network_mn, mip_optimizer)
    # return WM.optimize_model!(wm_master; relax_integrality = false)

    # control_settings = get_control_settings_from_result(result)
    return simulate_control_settings_sequential(wm_cq, control_settings)

    # sim_result = simulate!(network, heuristic_result, nlp_optimizer)
end


function solve_owf(network_path::String, obbt_optimizer, owf_optimizer, nlp_optimizer; use_new::Bool = true, kwargs...)
    network = WM.parse_file(network_path)
    network_mn = WM.make_multinetwork(network)
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)
    update_multinetwork_breakpoints!(network_mn, result_micp)

    # _update_tank_time_series_heur!(network, result_micp)
    wm_master = construct_owf_model(network_mn, owf_optimizer; use_new = use_new, kwargs...)
    # heuristic_setting = calc_heuristic(wm_master, network, obbt_optimizer, nlp_optimizer)
    # heur_cost = set_warm_start_from_setting!(wm_master, network, heuristic_setting, nlp_optimizer)
    # WM.Memento.info(LOGGER, "Heuristic found solution with objective $(round(heur_cost; digits = 2)).")

    # if use_new
    #     # Add binary-binary and binary-continuous pairwise cuts.
    #     WM.Memento.info(LOGGER, "Beginning cut preprocessing routine.")
    #     cut_time = @elapsed pairwise_cuts = compute_pairwise_cuts(network, obbt_optimizer)
    #     cut_time += @elapsed add_pairwise_cuts(wm_master, pairwise_cuts)
    #     cut_time += @elapsed add_pump_volume_cuts!(wm_master)
    #     WM.Memento.info(LOGGER, "Cut preprocessing completed in $(cut_time) seconds.")
    # end

    # # TODO: Remove this once Gurobi.jl interface is fixed.
    # wm_master.model.moi_backend.optimizer.model.has_generic_callback = false

    # # Add the lazy cut callback.
    # lazy_cut_stats = add_owf_lazy_cut_callback!(wm_master, network, heuristic_setting[1], nlp_optimizer)

    # if use_new
    #     # Add the user cut callback.
    #     WM.Memento.info(LOGGER, "Building and applying user cut callback routine.")
    #     user_cut_stats = add_owf_user_cut_callback!(wm_master)
    # end

    # Add the heuristic callback.
    # add_owf_heuristic_callback!(wm_master, network, nlp_optimizer, obbt_optimizer)

    # Optimize the master WaterModels model.
    WM.Memento.info(LOGGER, "Solving the master OWF problem.")
    # WM._relax_all_direction_variables!(wm_master)
    WM.JuMP.set_optimizer_attribute(wm_master.model, "TimeLimit", 0.0)
    result = WM.optimize_model!(wm_master; relax_integrality = false)

    set_branching_priorities!(wm_master)
    set_lazy_attributes!(wm_master)

    # add_price_cuts!(wm_master)

    WM.JuMP.set_optimizer_attribute(wm_master.model, "TimeLimit", 60.0)
    result = WM.optimize_model!(wm_master; relax_integrality = false)
    WM.Memento.info(LOGGER, "Solved for $(result["solve_time"]) seconds.")

    # # Print out relevant statistics of the solution process.
    # percent_lazy = round(lazy_cut_stats.time_elapsed / result["solve_time"] * 100.0; digits = 2)
    # WM.Memento.info(LOGGER, "Lazy cuts accounted for $(percent_lazy)% of solve time.")

    # if use_new
    #     percent_user = round(user_cut_stats.time_elapsed / result["solve_time"] * 100.0; digits = 2)
    #     WM.Memento.info(LOGGER, "User cuts accounted for $(percent_user)% of solve time.")
    # end

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
