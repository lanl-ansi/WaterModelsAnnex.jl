function solve_mn_owfh(file, model_constructor, optimizer; kwargs...)
    return WM.solve_model(file, model_constructor, optimizer, build_mn_owfh; multinetwork=true, kwargs...)
end


function run_mn_owfh(network, model_constructor, optimizer; kwargs...)
    Memento.warn(_LOGGER, "\"run_\" methods should be renamed \"solve_\" and will be deprecated in future versions.")
    return solve_mn_owfh(network, model_constructor, optimizer; kwargs...)
end

function build_mn_owfh(wm::WM.AbstractWaterModel)
    # Create head loss functions, if necessary.
    WM._function_head_loss(wm)

    for (n, network) in WM.nws(wm)
        # Physical variables.
        WM.variable_head(wm; nw=n)
        WM.variable_flow(wm; nw=n)
        WM.variable_pump_head_gain(wm; nw=n)
        WM.variable_pump_power(wm; nw=n)

        # Indicator (status) variables.
        WM.variable_des_pipe_indicator(wm; nw=n)
        WM.variable_pump_indicator(wm; nw=n)
        WM.variable_regulator_indicator(wm; nw=n)
        WM.variable_valve_indicator(wm; nw=n)

        # Create flow-related variables for node attachments.
        WM.variable_demand_flow(wm; nw=n)
        WM.variable_reservoir_flow(wm; nw=n)
        WM.variable_tank_flow(wm; nw=n)

        # Flow conservation at all nodes.
        for (i, node) in WM.ref(wm, :node; nw=n)
            WM.constraint_flow_conservation(wm, i; nw=n)
            WM.constraint_node_directionality(wm, i; nw=n)
        end

        # Constraints on pipe flows, heads, and physics.
        for (a, pipe) in WM.ref(wm, :pipe; nw=n)
            WM.constraint_pipe_flow(wm, a; nw=n)
            WM.constraint_pipe_head(wm, a; nw=n)
            WM.constraint_pipe_head_loss(wm, a; nw=n)
        end

        # Constraints on design pipe flows, heads, and physics.
        for (a, des_pipe) in WM.ref(wm, :des_pipe; nw=n)
            WM.constraint_on_off_des_pipe_flow(wm, a; nw=n)
            WM.constraint_on_off_des_pipe_head(wm, a; nw=n)
            WM.constraint_on_off_des_pipe_head_loss(wm, a; nw=n)
        end

        # Selection of design pipes along unique arcs.
        for (k, arc) in WM.ref(wm, :des_pipe_arc; nw=n)
            WM.constraint_des_pipe_flow(wm, k, arc[1], arc[2]; nw=n)
            WM.constraint_des_pipe_head(wm, k, arc[1], arc[2]; nw=n)
            WM.constraint_des_pipe_selection(wm, k, arc[1], arc[2]; nw=n)
        end

        # Constraints on pump flows, heads, and physics.
        for (a, pump) in WM.ref(wm, :pump; nw=n)
            WM.constraint_on_off_pump_head(wm, a; nw=n)
            WM.constraint_on_off_pump_head_gain(wm, a; nw=n)
            WM.constraint_on_off_pump_flow(wm, a; nw=n)
            WM.constraint_on_off_pump_power(wm, a; nw=n)
        end

        # Constraints on short pipe flows and heads.
        for (a, regulator) in WM.ref(wm, :regulator; nw=n)
            WM.constraint_on_off_regulator_head(wm, a; nw=n)
            WM.constraint_on_off_regulator_flow(wm, a; nw=n)
        end

        # Constraints on short pipe flows and heads.
        for (a, short_pipe) in WM.ref(wm, :short_pipe; nw=n)
            WM.constraint_short_pipe_head(wm, a; nw=n)
            WM.constraint_short_pipe_flow(wm, a; nw=n)
        end

        # Constraints on valve flows and heads.
        for (a, valve) in WM.ref(wm, :valve; nw=n)
            WM.constraint_on_off_valve_head(wm, a; nw=n)
            WM.constraint_on_off_valve_flow(wm, a; nw=n)
        end
    end

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(WM.nw_ids(wm)))

    # Start with the first network, representing the initial time step.
    n_1, n_f = network_ids[1], network_ids[end]

    # Constraints on tank volumes.
    for (i, tank) in WM.ref(wm, :tank; nw = n_1)
        # Set initial conditions of tanks.
        WM.constraint_tank_volume(wm, i; nw = n_1)
    end

    # # Constraints on tank volumes.
    # for n_2 in network_ids[2:end]
    #     # Constrain tank volumes after the initial time step.
    #     for (i, tank) in WM.ref(wm, :tank; nw = n_2)
    #         WM.constraint_tank_volume(wm, i, n_1, n_2)
    #     end

    #     n_1 = n_2 # Update the first network used for integration.
    # end

    for i in WM.ids(wm, n_f, :tank)
        WM.constraint_tank_volume_recovery(wm, i, n_1, n_f)
    end

    # Add the optimal water flow objective.
    WM.objective_owf(wm)
end