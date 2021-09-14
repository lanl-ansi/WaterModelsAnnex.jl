function run_obbt(inp_path::String, modification_path::String, cuts_path::String, time_limit::Float64, optimizer)
    # Read in the network and correct data.
    network = WM.parse_file(inp_path; skip_correct = true);
    modifications = WM.parse_file(modification_path; skip_correct = true);
    WM._IM.update_data!(network, modifications);
    WM.correct_network_data!(network);

    WM.set_flow_partitions_si!(network, 1.0, 1.0e-4);
    flow_partition_func = x -> WM.set_flow_partitions_si!(x, 1.0, 1.0e-4);

    solve_obbt_owf_switching!(network, optimizer; use_relaxed_network = true,
        model_type = WM.PWLRDWaterModel, time_limit = time_limit,
        cuts = load_pairwise_cuts(cuts_path), max_iter = 9999, solve_relaxed = false,
        flow_partition_func = flow_partition_func);

    return network
end