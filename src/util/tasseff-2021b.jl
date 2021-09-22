function run_obbt(inp_path::String, modification_path::String, cuts_path::String, time_limit::Float64, optimizer)
    # Read in the network and correct data.
    network = WM.parse_file(inp_path; skip_correct = true);
    modifications = WM.parse_file(modification_path; skip_correct = true);
    WM._IM.update_data!(network, modifications);
    WM.correct_network_data!(network);

    WM.set_flow_partitions_si!(network, 0.1, 1.0e-4);
    flow_partition_func = x -> WM.set_flow_partitions_si!(x, 0.1, 1.0e-4);

    solve_obbt_owf_switching!(network, optimizer; use_relaxed_network = true,
        model_type = WM.PWLRDWaterModel, time_limit = time_limit,
        cuts = load_pairwise_cuts(cuts_path), max_iter = 9999, solve_relaxed = false,
        flow_partition_func = flow_partition_func);

    return network
end


function run_obbt_mn_relaxed!(network_mn::Dict{String, <:Any}, cuts_path::String, time_limit::Float64, optimizer)
    # Set flow partitioning scheme.
    WM.set_flow_partitions_si!(network_mn, 100.0, 1.0e-4);
    flow_partition_func = x -> WM.set_flow_partitions_si!(x, 100.0, 1.0e-4);

    # Run a continuous relaxation-based OBBT on the data.
    solve_obbt_owf_switching!(network_mn, optimizer; use_relaxed_network = false,
        model_type = WM.LRDWaterModel, time_limit = time_limit,
        cuts = load_pairwise_cuts(cuts_path), max_iter = 9999, solve_relaxed = true,
        flow_partition_func = flow_partition_func);
end


function run_obbt_mn!(network_mn::Dict{String, <:Any}, cuts_path::String, time_limit::Float64, optimizer)
    # Set flow partitioning scheme.
    WM.set_flow_partitions_si!(network_mn, 100.0, 1.0e-4);
    flow_partition_func = x -> WM.set_flow_partitions_si!(x, 100.0, 1.0e-4);

    # Run a PWLRD-based OBBT on the data.
    solve_obbt_owf_switching!(network_mn, optimizer; use_relaxed_network = false,
        model_type = WM.LRDWaterModel, time_limit = time_limit,
        cuts = load_pairwise_cuts(cuts_path), max_iter = 9999, solve_relaxed = false,
        flow_partition_func = flow_partition_func);
end


function solve_owf_upper_bounds(network_mn::Dict, network::Dict, cuts_path::String, build_method::Function, formulation::Type, mip_optimizer, nlp_optimizer)
    # Solve a continuously-relaxed version of the problem.
    wm_micp = WM.instantiate_model(network_mn, formulation, build_method)
    WM.JuMP.set_optimizer(wm_micp.model, nlp_optimizer)
    WM.optimize_model!(wm_micp; relax_integrality = true)

    # Instantiate the model and add cutting planes.
    wm = WM.instantiate_model(network_mn, formulation, build_method)
    add_pairwise_cuts(wm, load_pairwise_cuts(cuts_path))
 
    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, mip_optimizer)
    #WM._MOI.set(wm.model, WM._MOI.NumberOfThreads(), 1)

    # TODO: Remove this once Gurobi.jl interface is fixed.
    #wm.model.moi_backend.optimizer.model.has_generic_callback = false

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
