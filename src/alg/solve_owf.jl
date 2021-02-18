function solve_obbt(network_path::String, modification_path::String, obbt_optimizer; time_limit::Float64 = 3600.0)
    # Read in the original network data.
    network = WM.parse_file(network_path; skip_correct = true)
    modifications = WM.parse_file(modification_path; skip_correct = true)
    WM._IM.update_data!(network, modifications)
    WM.correct_network_data!(network)

    # Tighten bounds of variables in the network.
    WM.solve_obbt_owf!(
        network, obbt_optimizer; model_type = WM.PWLRDWaterModel,
        solve_relaxed = false, time_limit = time_limit, max_iter = 1,
        ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10))
    
    # Get tightened network data.
    return network
end


function solve_owf(network_path::String, modification_path::String, obbt_optimizer, owf_optimizer, nlp_optimizer; kwargs...)
    # Tighten the bounds in the network.
    network = solve_obbt(network_path, modification_path, obbt_optimizer)
    result = solve_owf(network_path, network, obbt_optimizer, owf_optimizer, nlp_optimizer; kwargs...)
    return result
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


function construct_owf_model(network::Dict{String, Any}, owf_optimizer)
    # Construct the OWF model.
    network_mn = WM.make_multinetwork(network)
    ext = Dict(:pipe_breakpoints => 3, :pump_breakpoints => 3)
    #wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf; ext = ext)
    wm = WM.instantiate_model(network_mn, CDWaterModel, WM.build_mn_owf; ext = ext)

    # Constrain an auxiliary objective variable by the objective function.
    objective_function = WM.JuMP.objective_function(wm.model)
    objective_var = WM.JuMP.@variable(wm.model, base_name = "obj_aux", lower_bound = 0.0)
    WM.JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective_var)
    WM.JuMP.@constraint(wm.model, objective_function <= objective_var)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, owf_optimizer)
    #WM._MOI.set(wm.model, WM._MOI.NumberOfThreads(), 1)

    # Return the WaterModels object.
    return wm
end


function construct_relaxed_owf_model(network::Dict{String, Any}, nlp_optimizer)
    # Construct the relaxed OWF model.
    network_mn = WM.make_multinetwork(network)
    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf)
    WM.relax_all_binary_variables!(wm)

    # Set a feasibility-only objective.
    WM.JuMP.@objective(wm.model, WM._MOI.FEASIBILITY_SENSE, 0.0)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, nlp_optimizer)

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


function solve_owf(network_path::String, network, obbt_optimizer, owf_optimizer, nlp_optimizer; use_pairwise_cuts::Bool = true)
    wm = construct_owf_model(network, owf_optimizer)

    if use_pairwise_cuts
        pairwise_cuts = compute_pairwise_cuts(network, obbt_optimizer)
        add_pairwise_cuts(wm, pairwise_cuts)
    end

    # Solve the OWF optimization problem.
    result = WM.optimize_model!(wm)

    # Return the optimization result dictionary.
    return result
end
