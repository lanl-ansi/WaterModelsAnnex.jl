import LightGraphs

function generate_graph(wm::GenericWaterModel, q::Dict{Int, Float64}, h::Dict{Int, Float64}, n::Int=wm.cnw)
    graph = LightGraphs.SimpleDiGraph(length(wm.ref[:nw][n][:nodes]))
    node_map = sort(collect(ids(wm, n, :nodes)))

    for (i, node) in enumerate(node_map)
        arcs_fr = ref(wm, n, :node_arcs_fr, node)

        for (a, node_i, node_j) in arcs_fr
            j = findfirst(x -> x == node_j, node_map)

            if q[a] >= 0.0
                LightGraphs.add_edge!(graph, i, j)
            else
                LightGraphs.add_edge!(graph, j, i)
            end
        end
    end

    return graph
end

function add_upstream_cut!(wm::GenericWaterModel, graph::LightGraphs.SimpleDiGraph,
    q::Dict{Int, Float64}, resistance_indices::Dict{Int, Int}, i::Int, n::Int=wm.cnw)
    node_map = sort(collect(ids(wm, n, :nodes)))
    downstream_graph = LightGraphs.dfs_tree(graph, i)
    upstream_graph = LightGraphs.difference(graph, downstream_graph)
    lhs = JuMP.AffExpr(0.0)

    for edge in collect(LightGraphs.edges(upstream_graph))
        i, j = [node_map[edge.src], node_map[edge.dst]]
        a = WM.get_link_id(wm, i, j, n)
        r_id = resistance_indices[a]
        zero_indices = setdiff(1:length(wm.var[:nw][n][:x_res][a]), [r_id])
        lhs += wm.var[:nw][n][:x_res][a][r_id]
        lhs -= sum(wm.var[:nw][n][:x_res][a][zero_indices])
    end

    rhs = LightGraphs.ne(upstream_graph) - 1
    JuMP.@constraint(wm.model, lhs <= rhs)
end

function add_decomp_cut!(wm::GenericWaterModel, optimizer::JuMP.OptimizerFactory, n::Int=wm.cnw)
    add_no_good_cut!(wm, optimizer, n) # Always add a no-good cut.

    q, h, resistance_indices = get_cnlp_solution(wm, optimizer)
    qlb, qub, hlb, hub = check_solution_bounds(wm, q, h, resistance_indices)
    graph = generate_graph(wm, q, h, n)

    node_list = collect(keys(filter(d -> !d.second, hlb)))
    node_list = vcat(node_list, collect(keys(filter(d -> !d.second, hub))))

    for i in Set(node_list)
        add_upstream_cut!(wm, graph, q, resistance_indices, i, n)
    end

    arc_list = collect(keys(filter(d -> !d.second, qlb)))
    arc_list = vcat(arc_list, collect(keys(filter(d -> !d.second, qub))))

    for a in Set(arc_list)
        i = WM.ref(wm, n, :links, a)["node_fr"]
        j = WM.ref(wm, n, :links, a)["node_to"]

        if q[a] >= 0.0
            add_upstream_cut!(wm, graph, q, resistance_indices, i, n)
        else
            add_upstream_cut!(wm, graph, q, resistance_indices, j, n)
        end
    end
end
