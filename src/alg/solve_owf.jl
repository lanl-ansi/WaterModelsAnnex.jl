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
    # Tighten the bounds in the network.
    network = solve_obbt(network_path, modification_path, obbt_optimizer; kwargs...)
    return solve_owf(network_path, network, obbt_optimizer, owf_optimizer, nlp_optimizer; kwargs...)
end


function compute_pairwise_cuts(network::Dict{String, Any}, obbt_optimizer)
    # Relax the network.
    WM._relax_network!(network)

    # Get pairwise cutting planes from the network-relaxed problem.
    ext = Dict{Symbol, Any}(:pipe_breakpoints => 10, :pump_breakpoints => 10)
    wm = WM.instantiate_model(network, WM.PWLRDWaterModel, WM.build_owf; ext = ext)
    WM.JuMP.set_optimizer(wm.model, obbt_optimizer)
    problem_sets = WM._get_pairwise_problem_sets(wm)
    cuts = WM._compute_pairwise_cuts!(wm, problem_sets)

    # Unrelax the network.
    WM._fix_demands!(network)
    WM._fix_tanks!(network)
    WM._fix_reservoirs!(network)

    # Return the data structure comprising cuts.
    return cuts
end


function construct_owf_model(network::Dict{String, Any}, owf_optimizer; use_lrdx::Bool = true, kwargs...)
    # Construct the OWF model.
    network_mn = WM.make_multinetwork(network)
    ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10)
    wm = WM.instantiate_model(network_mn, WM.LRDWaterModel, WM.build_mn_owf; ext = ext)

    if use_lrdx
        wm = WM.instantiate_model(network_mn, LRDXWaterModel, WM.build_mn_owf; ext = ext)
    else
        wm = WM.instantiate_model(network_mn, WM.LRDWaterModel, WM.build_mn_owf; ext = ext)
    end

    # Constrain an auxiliary objective variable by the objective function.
    objective_function = WM.JuMP.objective_function(wm.model)
    objective_var = WM.JuMP.@variable(wm.model, base_name = "obj_aux", lower_bound = 0.0)
    WM.JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective_var)
    WM.JuMP.@constraint(wm.model, objective_function <= objective_var)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, owf_optimizer)
    WM._MOI.set(wm.model, WM._MOI.NumberOfThreads(), 1)

    # Return the WaterModels object.
    return wm
end


function add_pairwise_cuts(wm::WM.AbstractWaterModel, cuts::Array{WM._PairwiseCut, 1})
    # Add the pairwise cuts obtained from the relaxed problem to the OWF problem.
    for nw_id in WM.nw_ids(wm)
        # Use the same cuts for all subnetworks of the multinetwork.
        map(x -> x.variable_index_1.network_index = nw_id, cuts)
        map(x -> x.variable_index_2.network_index = nw_id, cuts)

        # Add the collection of pairwise cuts for the subnetwork.
        WM._add_pairwise_cuts!(wm, cuts)
    end
end


function solve_owf(network_path::String, network, obbt_optimizer, owf_optimizer, nlp_optimizer; use_pairwise_cuts::Bool = true, kwargs...)
    # Construct the OWF model that will serve as the master problem.
    wm = construct_owf_model(network, owf_optimizer; kwargs...)

    if use_pairwise_cuts
        # Add binary-binary and binary-continuous pairwise cuts.
        pairwise_cuts = compute_pairwise_cuts(network, obbt_optimizer)
        add_pairwise_cuts(wm, pairwise_cuts)
    end

    # Solve a relaxation of the master problem to begin.
    undo_relax = WM.JuMP.relax_integrality(wm.model)
    result_relaxed = WM.optimize_model!(wm)
    undo_relax() # Undo the continuous relaxation on the model.

    # TODO: Remove this once Gurobi.jl interface is fixed.
    wm.model.moi_backend.optimizer.model.has_generic_callback = false

    # Update the tank level time series to be used in finding an initial solution.
    _update_tank_time_series!(network, result_relaxed)

    # Find an initial feasible solution using a heuristic.
    initial_solution_time = @elapsed result_initial_solution =
        compute_initial_solution(network, obbt_optimizer, nlp_optimizer)

    # Warm start the primary WaterModels model.
    _set_initial_solution(wm, result_initial_solution)

    # Add the lazy cut callback.
    add_owf_lazy_cut_callback!(wm, nlp_optimizer)

    # Add the user cut callback.
    add_owf_user_cut_callback!(wm)

    # Add the heuristic callback.
    add_owf_heuristic_callback!(wm)

    # Optimize the master WaterModels model.
    return WM.optimize_model!(wm)
end


function _set_initial_solution(wm::WM.AbstractWaterModel, result::Dict{String, <:Any})
    for nw in WM.nw_ids(wm)
        for (i, reservoir) in WM.ref(wm, nw, :reservoir)
            qr = WM.var(wm, nw, :q_reservoir, i)
            qr_sol = result["solution"]["nw"][string(nw)]["reservoir"][string(i)]["q"]
            JuMP.set_start_value(qr, qr_sol)
        end

        for (i, tank) in WM.ref(wm, nw, :tank)
            qt = WM.var(wm, nw, :q_tank, i)
            qt_sol = result["solution"]["nw"][string(nw)]["tank"][string(i)]["q"]
            JuMP.set_start_value(qt, qt_sol)
        end

        for (i, pump) in WM.ref(wm, nw, :pump)
            Ps = WM.var(wm, nw, :Ps_pump, i)
            Ps_sol = result["solution"]["nw"][string(nw)]["pump"][string(i)]["Ps"]
            JuMP.set_start_value(Ps, Ps_sol)

            g = WM.var(wm, nw, :g_pump, i)
            g_sol = result["solution"]["nw"][string(nw)]["pump"][string(i)]["g"]
            JuMP.set_start_value(g, g_sol)
        end

        for comp_type in [:pipe, :des_pipe]
            for (i, comp) in WM.ref(wm, nw, comp_type)
                comp_sol = result["solution"]["nw"][string(nw)][string(comp_type)][string(i)]
                dhp = WM.var(wm, nw, Symbol("dhp_" * string(comp_type)), i)
                dhn = WM.var(wm, nw, Symbol("dhn_" * string(comp_type)), i)

                h_i = result["solution"]["nw"][string(nw)]["node"][string(comp["node_fr"])]["h"]
                h_j = result["solution"]["nw"][string(nw)]["node"][string(comp["node_to"])]["h"]

                dhp_sol, dhn_sol = max(0.0, h_i - h_j), max(0.0, h_j - h_i)
                JuMP.set_start_value(dhp, dhp_sol)
                JuMP.set_start_value(dhn, dhn_sol)
            end
        end
        
        for (i, node) in WM.ref(wm, nw, :node)
            h = WM.var(wm, nw, :h, i)
            h_sol = result["solution"]["nw"][string(nw)]["node"][string(i)]["h"]
            JuMP.set_start_value(h, h_sol)
        end

        for comp_type in [:pipe, :pump, :regulator, :short_pipe, :valve]
            for (i, comp) in WM.ref(wm, nw, comp_type)
                y = WM.var(wm, nw, Symbol("y_" * string(comp_type)), i)
                qp = WM.var(wm, nw, Symbol("qp_" * string(comp_type)), i)
                qn = WM.var(wm, nw, Symbol("qn_" * string(comp_type)), i)

                comp_sol = result["solution"]["nw"][string(nw)][string(comp_type)][string(i)]
                JuMP.set_start_value(qp, max(0.0, comp_sol["q"]))
                JuMP.set_start_value(qn, max(0.0, -comp_sol["q"]))
                y_start = comp_sol["q"] > 1.0e-4 ? 1.0 : 0.0
                JuMP.set_start_value(y, y_start)
            end
        end

        for comp_type in [:pump, :regulator, :valve]
            for (i, comp) in WM.ref(wm, nw, comp_type)
                z = WM.var(wm, nw, Symbol("z_" * string(comp_type)), i)
                comp_sol = result["solution"]["nw"][string(nw)][string(comp_type)][string(i)]
                z_start = Float64(Int(comp_sol["status"]))
                JuMP.set_start_value(z, z_start)
            end
        end
    end
end


function compute_initial_solution(network::Dict{String, <:Any}, obbt_optimizer, nlp_optimizer)
    schedules = calc_possible_schedules(network, obbt_optimizer, nlp_optimizer)
    weights = solve_heuristic_problem(network, schedules, obbt_optimizer)
    return solve_heuristic_master(network, schedules, weights, nlp_optimizer, obbt_optimizer)
end