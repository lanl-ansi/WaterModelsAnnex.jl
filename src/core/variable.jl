function variable_volume_sum(wm::WM.AbstractWaterModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables for total hydraulic head.
    V_sum = WM.var(wm, nw)[:V_sum] = JuMP.@variable(wm.model,
        [i in WM.ids(wm, nw, :tank_group)], base_name = "$(nw)_V_sum",
        start = WM.comp_start_value(WM.ref(wm, nw, :tank_group, i), "V_sum_start"))
 
    if bounded
        for (i, tank_group) in WM.ref(wm, nw, :tank_group)
            if haskey(tank_group, "V_sum_min")
                V_sum_min = max(0.0, tank_group["V_sum_min"])
            else
                V_sum_min = 0.0

                for tank in WM.ref.(Ref(wm), nw, :tank, tank_group["tank_indices"]) 
                    node = WM.ref(wm, nw, :node, tank["node"])
                    level_min = node["head_min"] - node["elevation"]
                    V_sum_min += 0.25 * pi * level_min * tank["diameter"]^2
                end
            end

            if haskey(tank_group, "V_sum_max")
                V_sum_max = max(0.0, tank_group["V_sum_max"])
            else
                V_sum_max = 0.0

                for tank in WM.ref.(Ref(wm), nw, :tank, tank_group["tank_indices"]) 
                    node = WM.ref(wm, nw, :node, tank["node"])
                    level_max = node["head_max"] - node["elevation"]
                    V_sum_max += 0.25 * pi * level_max * tank["diameter"]^2
                end
            end

            # Set the lower and upper bounds for each head.
            JuMP.set_lower_bound(V_sum[i], V_sum_min)
            JuMP.set_upper_bound(V_sum[i], V_sum_max)

            # Set start value for the head variable with possibly better data.
            V_sum_mid = 0.5 * (V_sum_min + V_sum_max)
            V_sum_start = WM.comp_start_value(tank_group, "V_sum_start", V_sum_mid)
            JuMP.set_start_value(V_sum[i], V_sum_start)
        end
    end

    report && WM.sol_component_value(wm, nw, :tank_group,
        :V_sum, WM.ids(wm, nw, :tank_group), V_sum)
end
