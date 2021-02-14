function WM.build_wf(wm::AbstractCDModel)
    # Create head loss functions, if necessary.
    WM._function_head_loss(wm)

    # Physical variables.
    variable_flow(wm; bounded = false)

    # Create flow-related variables for node attachments.
    WM.variable_demand_flow(wm)
    WM.variable_reservoir_flow(wm)
    WM.variable_tank_flow(wm)

    # Flow conservation at all nodes.
    for (i, node) in WM.ref(wm, :node)
        WM.constraint_flow_conservation(wm, i)
    end

    # Add the objective.
    WM.objective_wf(wm)
end