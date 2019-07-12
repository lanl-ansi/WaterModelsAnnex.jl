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

    wf = WMs.build_generic_model(network, WMs.CNLPWaterModel, WMs.get_post_wf(alpha))
    wf_result = WMs.solve_generic_model(wf, nlp_optimizer)

    q = Dict(parse(Int, a) => pipe["q"] for (a, pipe) in wf_result["solution"]["pipes"])
    r = Dict(parse(Int, a) => pipe["r"] for (a, pipe) in wf_result["solution"]["pipes"])

    junction_ids = collect(ids(wf, n, :junctions))
    reservoir_ids = collect(ids(wf, n, :reservoirs))
    link_ids = collect(ids(wf, n, :links))
    node_ids = [junction_ids; reservoir_ids]
    node_mapping = Dict{Int, Int}(node_ids[i] => i for i in 1:length(node_ids))

    num_reservoirs = length(reservoir_ids)
    num_nodes = length(node_ids)
    num_links = length(link_ids)

    # Create matrices for the left- and right-hand sides (for Ax = b).
    A = zeros(Float64, num_links + num_reservoirs, num_nodes)
    b = zeros(Float64, num_links + num_reservoirs, 1)

    for (row, a) in enumerate(link_ids)
        link = wf.ref[:nw][n][:links][a]

        node_i = link["f_id"]
        A[row, node_mapping[node_i]] = 1.0

        node_j = link["t_id"]
        A[row, node_mapping[node_j]] = -1.0

        L = link["length"]
        b[row] = L * r[a] * q[a] * abs(q[a])^(alpha - 1.0)
    end

    for (i, reservoir_id) in enumerate(reservoir_ids)
        A[num_links + i, node_mapping[reservoir_id]] = 1
        b[num_links + i] = wf.ref[:nw][n][:reservoirs][reservoir_id]["head"]
    end

    h = A \ b # Get the solution for head variables.
    h = Dict{Int, Float64}(i => h[node_mapping[i]] for i in node_ids)

    return q, h, resistance_indices
end
