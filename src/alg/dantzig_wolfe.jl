function solve_dantzig_wolfe(network_path::String, modification_path::String, optimizer)
    coluna = JuMP.optimizer_with_attributes(Coluna.Optimizer,
        "params" => Coluna.Params(solver =
        Coluna.Algorithm.TreeSearchAlgorithm()),
        "default_optimizer" => optimizer)
    
    # Set up the Block Decomposition model.
    model = BD.BlockModel(coluna, automatic_decomposition = true)
    # BD.objectiveprimalbound!(model, 1188.9444282)
    #Coluna.set_initial_primal_bound!(model, 1188.9444282);
    #model = JuMP.Model(optimizer);

    #BD.@dantzig_wolfe_decomposition(model, decomposition)

    # Initialize the network and model.
    network = WM.parse_file(network_path)
    network_mn = WM.make_multinetwork(network)
    ext = Dict(:pipe_breakpoints => 5, :pump_breakpoints => 5)
    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel,
        WM.build_mn_owf; ext = ext, jump_model = model)
    
    return wm

    # Solve the optimization problem.
    JuMP.optimize!(wm.model)
end