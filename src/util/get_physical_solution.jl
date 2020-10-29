function get_cnlp_network(wm::GenericWaterModel, n::Int=wm.cnw)
    # Create a copy of the original network data.
    network = deepcopy(wm.data)

    for (a, link) in wm.ref[:nw][n][:links_ne]
        # Remove unnecessary data from the network specification.
        delete!(network["pipes"][string(a)], "diameter")
        delete!(network["pipes"][string(a)], "diameters")
        delete!(network["pipes"][string(a)], "resistances")
    end

    # Return the network.
    return network
end

function get_cnlp_solution(wm::GenericWaterModel, nlp_optimizer::JuMP.OptimizerFactory, n::Int=wm.cnw)
    alpha = wm.ref[:nw][n][:alpha]
    network = get_cnlp_network(wm, n)
    resistance_indices = Dict{Int, Int}(a => 0 for a in collect(ids(wm, n, :links)))

    for (a, link) in wm.ref[:nw][n][:links_ne]
        x_res, r_id = findmax(JuMP.value.(wm.var[:nw][wm.cnw][:x_res][a]))
        resistance_indices[a] = r_id
        resistance = wm.ref[:nw][n][:resistance][a][r_id]
        network["pipes"][string(a)]["resistance"] = resistance
    end

    wf = WM.build_generic_model(network, WM.CNLPWaterModel, WM.post_wf)
    wf_result = WM.solve_generic_model(wf, nlp_optimizer)

    q = Dict(parse(Int, a) => pipe["q"] for (a, pipe) in wf_result["solution"]["pipes"])
    r = Dict(parse(Int, a) => pipe["r"] for (a, pipe) in wf_result["solution"]["pipes"])
    h = Dict(parse(Int, i) => node["h"] for (i, node) in wf_result["solution"]["nodes"])

    return q, h, resistance_indices, wf_result["objective_value"]
end
