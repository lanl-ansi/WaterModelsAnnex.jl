function get_cnlp_network(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Create a copy of the original network data.
    network = deepcopy(wm.data)

    for (a, link) in wm.ref[:nw][n][:link_ne]
        # Remove unnecessary data from the network specification.
        delete!(network["pipe"][string(a)], "diameter")
        delete!(network["pipe"][string(a)], "diameters")
        delete!(network["pipe"][string(a)], "resistances")
    end

    # Return the network.
    return network
end

function get_cnlp_solution(wm::AbstractWaterModel, nlp_optimizer::JuMP.OptimizerFactory, n::Int=wm.cnw)
    alpha = wm.ref[:nw][n][:alpha]
    network = get_cnlp_network(wm, n)
    resistance_indices = Dict{Int, Int}(a => 0 for a in collect(ids(wm, n, :link)))

    for (a, link) in wm.ref[:nw][n][:link_ne]
        x_res, r_id = findmax(JuMP.value.(wm.var[:nw][n][:x_res][a]))
        resistance_indices[a] = r_id
        resistance = wm.ref[:nw][n][:resistance][a][r_id]
        network["pipe"][string(a)]["resistance"] = resistance
    end

    wf = WMs.build_model(network, WMs.CNLPWaterModel, WMs.post_wf)
    wf_result = WMs.optimize_model!(wf, nlp_optimizer)

    q = Dict(parse(Int, a) => pipe["q"] for (a, pipe) in wf_result["solution"]["pipe"])
    r = Dict(parse(Int, a) => pipe["r"] for (a, pipe) in wf_result["solution"]["pipe"])
    h = Dict(parse(Int, i) => node["h"] for (i, node) in wf_result["solution"]["node"])

    return q, h, resistance_indices, wf_result["objective"]
end
