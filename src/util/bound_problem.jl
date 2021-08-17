function _get_bound_problems_tank_groups(wm::WM.AbstractWaterModel; limit::Bool = false)
    return vcat(_get_bound_problems_tank_groups.(Ref(wm), WM.nw_ids(wm); limit = limit)...)
end


function _get_bound_problems_tank_groups(wm::WM.AbstractWaterModel, nw::Int; limit::Bool = false)
    return vcat(_get_bound_problems_tank_group.(Ref(wm),
        WM.ids(wm, nw, :tank_group), nw; limit = limit)...)
end


function _get_bound_problems_tank_group(wm::WM.AbstractWaterModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(WM.var(wm, nw), :V_sum) && i in [x[1] for x in eachindex(WM.var(wm, nw, :V_sum))]
        V_sum_vid = WM._VariableIndex(nw, :tank_group, :V_sum, i)

        V_sum_min = WM._get_lower_bound_from_index(wm, V_sum_vid)
        bp_V_sum_min = WM.BoundProblem(WM._MOI.MIN_SENSE, V_sum_vid, [],
            [], "V_sum_min", V_sum_min, 1.0e-3, true)

        V_sum_max = WM._get_upper_bound_from_index(wm, V_sum_vid)
        bp_V_sum_max = WM.BoundProblem(WM._MOI.MAX_SENSE, V_sum_vid, [],
            [], "V_sum_max", V_sum_max, 1.0e-3, true)

        if limit
            return Vector{WM.BoundProblem}([])
        else
            return Vector{WM.BoundProblem}([bp_V_sum_min, bp_V_sum_max])
        end
    else
        return Vector{WM.BoundProblem}([])
    end
end


function _get_bound_problems(wm::WM.AbstractWaterModel; limit::Bool = false)::Vector{WM.BoundProblem}
    # Create the sets of bound-tightening problems.
    bps_node = WM._get_bound_problems_nodes(wm; limit = limit)
    bps_pipe = WM._get_bound_problems_pipes(wm; limit = limit)
    bps_pump = WM._get_bound_problems_pumps(wm; limit = limit)
    bps_regulator = WM._get_bound_problems_regulators(wm; limit = limit)
    bps_short_pipe = WM._get_bound_problems_short_pipes(wm; limit = limit)
    bps_valve = WM._get_bound_problems_valves(wm; limit = limit)
    bps_tank_groups = _get_bound_problems_tank_groups(wm; limit = limit)
    return vcat(bps_pump, bps_valve, bps_regulator,
        bps_pipe, bps_short_pipe, bps_node, bps_tank_groups)
end