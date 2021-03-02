function WM.build_wf(wm::AbstractCDModel)
    # Create head loss functions, if necessary.
    WM._function_head_loss(wm)

    # Physical variables.
    WM.variable_flow(wm)
    WM.variable_head(wm)
    WM.variable_pump_head_gain(wm)
    WM.variable_pump_power(wm)

    # Indicator (status) variables.
    WM.variable_des_pipe_indicator(wm)
    WM.variable_pump_indicator(wm)
    WM.variable_regulator_indicator(wm)
    WM.variable_valve_indicator(wm)

    # Create flow-related variables for node attachments.
    WM.variable_demand_flow(wm)
    WM.variable_reservoir_flow(wm)
    WM.variable_tank_flow(wm)

    # Flow conservation at all nodes.
    for (i, node) in WM.ref(wm, :node)
        WM.constraint_flow_conservation(wm, i)
    end

    # Constraints on pipe flows, heads, and physics.
    for (a, pipe) in WM.ref(wm, :pipe)
        WM.constraint_pipe_head(wm, a)
    end

    # Constraints on design pipe flows, heads, and physics.
    for (a, des_pipe) in WM.ref(wm, :des_pipe)
        WM.constraint_on_off_des_pipe_head(wm, a)
        WM.constraint_on_off_des_pipe_flow(wm, a)
    end

    # Selection of design pipes along unique arcs.
    for (k, arc) in WM.ref(wm, :des_pipe_arc)
        WM.constraint_des_pipe_head(wm, k, arc[1], arc[2])
        WM.constraint_des_pipe_selection(wm, k, arc[1], arc[2])
    end

    # Constraints on pump flows and heads.
    for (a, pump) in WM.ref(wm, :pump)
        WM.constraint_on_off_pump_head(wm, a)
        WM.constraint_on_off_pump_flow(wm, a)
        WM.constraint_on_off_pump_power(wm, a)
    end

    # Constraints on short pipe flows and heads.
    for (a, regulator) in WM.ref(wm, :regulator)
        WM.constraint_on_off_regulator_head(wm, a)
        WM.constraint_on_off_regulator_flow(wm, a)
    end

    # Constraints on short pipe flows and heads.
    for (a, short_pipe) in WM.ref(wm, :short_pipe)
        WM.constraint_short_pipe_head(wm, a)
        WM.constraint_short_pipe_flow(wm, a)
    end

    # Constraints on valve flows and heads.
    for (a, valve) in WM.ref(wm, :valve)
        WM.constraint_on_off_valve_head(wm, a)
        WM.constraint_on_off_valve_flow(wm, a)
    end

    # Constraints on tank volumes.
    for (i, tank) in WM.ref(wm, :tank)
        # Set the initial tank volume.
        WM.constraint_tank_volume(wm, i)
    end

    # Add the strong duality constraint.
    constraint_strong_duality(wm)

    # Add the objective.
    # WM.objective_wf(wm)
end


function WM.build_mn_wf(wm::AbstractCDModel)
    # Create head loss functions, if necessary.
    WM._function_head_loss(wm)

    for (n, network) in WM.nws(wm)
        # Physical variables.
        WM.variable_flow(wm; nw = n, bounded = true)
        WM.variable_head(wm; nw = n, bounded = true)
        WM.variable_pump_head_gain(wm; nw = n, bounded = true)
        WM.variable_pump_power(wm; nw = n, bounded = true)

        # Indicator (status) variables.
        WM.variable_des_pipe_indicator(wm; nw = n)
        WM.variable_pump_indicator(wm; nw = n)
        WM.variable_regulator_indicator(wm; nw = n)
        WM.variable_valve_indicator(wm; nw = n)

        # Create flow-related variables for node attachments.
        WM.variable_demand_flow(wm; nw = n)
        WM.variable_reservoir_flow(wm; nw = n)
        WM.variable_tank_flow(wm; nw = n)

        # Flow conservation at all nodes.
        for (i, node) in WM.ref(wm, n, :node)
            WM.constraint_flow_conservation(wm, i; nw = n)
        end

        # Constraints on pipe flows, heads, and physics.
        for (a, pipe) in WM.ref(wm, n, :pipe)
            WM.constraint_pipe_head(wm, a; nw = n)
        end

        # Constraints on design pipe flows and heads.
        for (a, des_pipe) in WM.ref(wm, n, :des_pipe)
            WM.constraint_on_off_des_pipe_head(wm, a; nw = n)
            WM.constraint_on_off_des_pipe_flow(wm, a; nw = n)
        end

        # Selection of design pipes along unique arcs.
        for (k, arc) in WM.ref(wm, n, :des_pipe_arc)
            WM.constraint_des_pipe_head(wm, k, arc[1], arc[2]; nw = n)
            WM.constraint_des_pipe_selection(wm, k, arc[1], arc[2]; nw = n)
        end

        # Constraints on pump flows and heads.
        for (a, pump) in WM.ref(wm, n, :pump)
            WM.constraint_on_off_pump_head(wm, a; nw = n)
            WM.constraint_on_off_pump_flow(wm, a; nw = n)
            WM.constraint_on_off_pump_power(wm, a; nw = n)
        end

        # Constraints on short pipe flows and heads.
        for (a, regulator) in WM.ref(wm, n, :regulator)
            WM.constraint_on_off_regulator_head(wm, a; nw = n)
            WM.constraint_on_off_regulator_flow(wm, a; nw = n)
        end

        # Constraints on short pipe flows and heads.
        for (a, short_pipe) in WM.ref(wm, n, :short_pipe)
            WM.constraint_short_pipe_head(wm, a; nw = n)
            WM.constraint_short_pipe_flow(wm, a; nw = n)
        end

        # Constraints on valve flows and heads.
        for (a, valve) in WM.ref(wm, n, :valve)
            WM.constraint_on_off_valve_head(wm, a; nw = n)
            WM.constraint_on_off_valve_flow(wm, a; nw = n)
        end
    end

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(WM.nw_ids(wm)))

    # Start with the first network, representing the initial time step.
    n_1 = network_ids[1]

    # Constraints on tank volumes.
    for (i, tank) in WM.ref(wm, :tank; nw = n_1)
        # Set initial conditions of tanks.
        WM.constraint_tank_volume(wm, i; nw = n_1)
    end

    # Constraints on tank volumes.
    for n_2 in network_ids[2:end]
        # Constrain tank volumes after the initial time step.
        for (i, tank) in WM.ref(wm, :tank; nw = n_2)
            WM.constraint_tank_volume(wm, i, n_1, n_2)
        end

        n_1 = n_2 # Update the first network used for integration.
    end

    # Add the strong duality constraint.
    constraint_strong_duality(wm)

    # Add the objective.
    #WM.objective_wf(wm)
end


function WM.build_wf(wm::Union{AbstractCDXModel, AbstractLRDXModel})
    # Create head loss functions, if necessary.
    WM._function_head_loss(wm)

    # Physical variables.
    WM.variable_head(wm)
    WM.variable_flow(wm)
    WM.variable_pump_head_gain(wm)
    WM.variable_pump_power(wm)
    variable_pipe_flow_nonlinear(wm)

    # Indicator (status) variables.
    WM.variable_des_pipe_indicator(wm)
    WM.variable_pump_indicator(wm)
    WM.variable_regulator_indicator(wm)
    WM.variable_valve_indicator(wm)

    # Create flow-related variables for node attachments.
    WM.variable_demand_flow(wm)
    WM.variable_reservoir_flow(wm)
    WM.variable_tank_flow(wm)

    # Flow conservation at all nodes.
    for (i, node) in WM.ref(wm, :node)
        WM.constraint_flow_conservation(wm, i)
        WM.constraint_node_directionality(wm, i)
    end

    # Constraints on pipe flows, heads, and physics.
    for (a, pipe) in WM.ref(wm, :pipe)
        WM.constraint_pipe_head(wm, a)
        WM.constraint_pipe_head_loss(wm, a)
        WM.constraint_pipe_flow(wm, a)
        constraint_pipe_head_loss_integrated(wm, a)
    end

    # Selection of design pipes along unique arcs.
    for (k, arc) in WM.ref(wm, :des_pipe_arc)
        WM.constraint_des_pipe_flow(wm, k, arc[1], arc[2])
        WM.constraint_des_pipe_head(wm, k, arc[1], arc[2])
        WM.constraint_des_pipe_selection(wm, k, arc[1], arc[2])
    end

    # Constraints on design pipe flows, heads, and physics.
    for (a, des_pipe) in WM.ref(wm, :des_pipe)
        WM.constraint_on_off_des_pipe_head(wm, a)
        WM.constraint_on_off_des_pipe_head_loss(wm, a)
        WM.constraint_on_off_des_pipe_flow(wm, a)
    end

    # Constraints on pump flows, heads, and physics.
    for (a, pump) in WM.ref(wm, :pump)
        WM.constraint_on_off_pump_head(wm, a)
        WM.constraint_on_off_pump_head_gain(wm, a)
        WM.constraint_on_off_pump_flow(wm, a)
        WM.constraint_on_off_pump_power(wm, a)
    end

    # Constraints on short pipe flows and heads.
    for (a, regulator) in WM.ref(wm, :regulator)
        WM.constraint_on_off_regulator_head(wm, a)
        WM.constraint_on_off_regulator_flow(wm, a)
    end

    # Constraints on short pipe flows and heads.
    for (a, short_pipe) in WM.ref(wm, :short_pipe)
        WM.constraint_short_pipe_head(wm, a)
        WM.constraint_short_pipe_flow(wm, a)
    end

    # Constraints on tank volumes.
    for (i, tank) in WM.ref(wm, :tank)
        # Set the initial tank volume.
        WM.constraint_tank_volume(wm, i)
    end

    # Constraints on valve flows and heads.
    for (a, valve) in WM.ref(wm, :valve)
        WM.constraint_on_off_valve_head(wm, a)
        WM.constraint_on_off_valve_flow(wm, a)
    end

    # Add the strong duality constraint.
    WM.constraint_strong_duality(wm)

    # Add the objective.
    WM.objective_wf(wm)
end


function WM.build_mn_wf(wm::Union{AbstractCDXModel, AbstractLRDXModel})
    # Create head loss functions, if necessary.
    WM._function_head_loss(wm)

    for (n, network) in WM.nws(wm)
        # Physical variables.
        WM.variable_head(wm; nw=n)
        WM.variable_flow(wm; nw=n)
        WM.variable_pump_head_gain(wm; nw=n)
        WM.variable_pump_power(wm; nw=n)

        # Additional variables for nonlinearities.
        variable_pipe_flow_nonlinear(wm; nw=n)
        variable_pipe_head_difference_nonlinear(wm; nw=n)

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

            constraint_pipe_flow_nonlinear(wm, a; nw=n)
            constraint_pipe_head_nonlinear(wm, a; nw=n)
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
    n_1 = network_ids[1]

    # Constraints on tank volumes.
    for (i, tank) in WM.ref(wm, :tank; nw = n_1)
        # Set initial conditions of tanks.
        WM.constraint_tank_volume(wm, i; nw = n_1)
    end

    # Constraints on tank volumes.
    for n_2 in network_ids[2:end]
        # Constrain tank volumes after the initial time step.
        for (i, tank) in WM.ref(wm, :tank; nw = n_2)
            WM.constraint_tank_volume(wm, i, n_1, n_2)
        end

        n_1 = n_2 # Update the first network used for integration.
    end

    # Add the strong duality constraint.
    constraint_strong_duality(wm)

    # Add the objective.
    WM.objective_wf(wm)
end