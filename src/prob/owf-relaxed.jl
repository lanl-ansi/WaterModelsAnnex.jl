function build_mn_owf_excess(wm::WM.AbstractWaterModel)
    # Build the water flow problem.
    build_mn_wf_excess(wm)

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(WM.nw_ids(wm)))

    # Ensure tanks recover their initial volume.
    n_1, n_f = network_ids[1], network_ids[end]

    for i in WM.ids(wm, n_f, :tank)
        WM.constraint_tank_volume_recovery(wm, i, n_1, n_f)
    end

    # Add the optimal water flow objective.
    WM.objective_owf(wm)
end
