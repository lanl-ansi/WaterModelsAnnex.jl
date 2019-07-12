function check_solution_bounds(wm::GenericWaterModel,
                               q::Dict{Int, Float64},
                               h::Dict{Int, Float64},
                               resistance_indices::Dict{Int, Int},
                               n::Int = wm.cnw)
    # Initialize dictionaries used to store arc infeasibility results.
    link_ids = collect(ids(wm, n, :links))
    q_sat_lb = Dict{Int, Bool}(a => true for a in link_ids)
    q_sat_ub = Dict{Int, Bool}(a => true for a in link_ids)

    # Compute bound satisfaction results for flow variables.
    for (a, link) in wm.ref[:nw][n][:links]
        # Get the selected resistance index for this arc.
        r_a = resistance_indices[a]

        if q[a] >= 0.0
            # Compute bound satisfaction for flow from i to j.
            q_sat_lb[a] = q[a] >= JuMP.lower_bound(wm.var[:nw][n][:qp_ne][a][r_a]) - 1.0e-7
            q_sat_ub[a] = q[a] <= JuMP.upper_bound(wm.var[:nw][n][:qp_ne][a][r_a]) + 1.0e-7
        else
            # Compute bound satisfaction for flow from j to i.
            q_sat_lb[a] = -q[a] >= JuMP.lower_bound(wm.var[:nw][n][:qn_ne][a][r_a]) - 1.0e-7
            q_sat_ub[a] = -q[a] <= JuMP.upper_bound(wm.var[:nw][n][:qn_ne][a][r_a]) + 1.0e-7
        end
    end

    # Initialize dictionaries used to store node infeasibility results.
    junction_ids = collect(ids(wm, n, :junctions))
    reservoir_ids = collect(ids(wm, n, :reservoirs))
    node_ids = [junction_ids; reservoir_ids]
    h_sat_lb = Dict{Int, Bool}(i => true for i in node_ids)
    h_sat_ub = Dict{Int, Bool}(i => true for i in node_ids)

    # Compute bound satisfaction results for head variables.
    for (i, junction) in wm.ref[:nw][n][:junctions]
        h_sat_lb[i] = h[i] >= JuMP.lower_bound(wm.var[:nw][n][:h][i]) - 1.0e-7
        h_sat_ub[i] = h[i] <= JuMP.upper_bound(wm.var[:nw][n][:h][i]) + 1.0e-7
    end

    # Return dictionaries of variable bound satisfaction results.
    return q_sat_lb, q_sat_ub, h_sat_lb, h_sat_ub
end
