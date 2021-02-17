function WM.build_wf(wm::AbstractCDModel)
    # Create head loss functions, if necessary.
    _function_head_loss(wm)

    # Physical variables.
    variable_flow(wm; bounded = true)
    WM.variable_head(wm; bounded = true)
    WM.variable_pump_head_gain(wm; bounded = true)
    WM.variable_pump_power(wm; bounded = true)

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

    # Selection of design pipes along unique arcs.
    for (k, arc) in WM.ref(wm, :des_pipe_arc)
        WM.constraint_des_pipe_head(wm, k, arc[1], arc[2])
        WM.constraint_des_pipe_selection(wm, k, arc[1], arc[2])
    end

    # Constraints on design pipe flows, heads, and physics.
    for (a, des_pipe) in WM.ref(wm, :des_pipe)
        WM.constraint_on_off_des_pipe_head(wm, a)
        WM.constraint_on_off_des_pipe_flow(wm, a)
    end

    # Constraints on pump flows and heads.
    for (a, pump) in WM.ref(wm, :pump)
        WM.constraint_on_off_pump_head(wm, a)
        WM.constraint_on_off_pump_flow(wm, a)
    end

    # Add the objective.
    WM.objective_wf(wm)
end