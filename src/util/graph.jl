function create_graph(data::Dict{String,<:Any})
    graph = MetaGraphs.MetaDiGraph(length(data["node"]))
    node_names = sort(collect(keys(data["node"])))
    node_map = Dict{String,Int}(x => i for (i, x) in enumerate(node_names))

    for node in values(data["node"])
        index = node_map[string(node["index"])]
        node_sym = Dict(Symbol(key) => value for (key, value) in node)
        MetaGraphs.set_props!(graph, index, node_sym)

        MetaGraphs.set_prop!(graph, index, :has_demand, false)
        MetaGraphs.set_prop!(graph, index, :has_reservoir, false)
        MetaGraphs.set_prop!(graph, index, :has_tank, false)
    end

    for demand in values(data["demand"])
        index = node_map[string(demand["node"])]
        MetaGraphs.set_prop!(graph, index, :has_demand, true) 
        demand_sym = Dict(Symbol(key) => value for (key, value) in demand)
        props = merge(demand_sym, MetaGraphs.props(graph, index))
        MetaGraphs.set_props!(graph, index, props)
    end

    for reservoir in values(data["reservoir"])
        index = node_map[string(reservoir["node"])]
        MetaGraphs.set_prop!(graph, index, :has_reservoir, true)
        reservoir_sym = Dict(Symbol(key) => value for (key, value) in reservoir)
        props = merge(reservoir_sym, MetaGraphs.props(graph, index))
        MetaGraphs.set_props!(graph, index, props)
    end

    for tank in values(data["tank"])
        index = node_map[string(tank["node"])]
        MetaGraphs.set_prop!(graph, index, :has_tank, true)
        tank_sym = Dict(Symbol(key) => value for (key, value) in tank)
        props = merge(tank_sym, MetaGraphs.props(graph, index))
        MetaGraphs.set_props!(graph, index, props)
    end

    for component_type in ["pipe", "pump", "regulator", "short_pipe", "valve"]
        for comp in values(data[component_type])
            i, j = node_map[string(comp["node_fr"])], node_map[string(comp["node_to"])]
            comp_sym = Dict(Symbol(key) => value for (key, value) in comp)
            comp_sym[:component_type] = component_type

            if comp["flow_max"] > 0.0
                LightGraphs.add_edge!(graph, i, j)
                MetaGraphs.set_props!(graph, i, j, comp_sym)
            end

            if comp["flow_min"] < 0.0
                LightGraphs.add_edge!(graph, j, i)
                MetaGraphs.set_props!(graph, j, i, comp_sym)
            end
        end
    end

    return (graph, node_map)
end


function update_graph_from_simulation_result!(wm::AbstractCQModel, control_setting::ControlSetting, graph)
    graph_nw, node_map = graph[1], graph[2]
    wm_data = WM.get_wm_data(wm.data)

    # for (k, vid) in enumerate(control_setting.variable_indices)
    #     val = control_setting.vals[k]
    #     comp = WM.ref(wm, vid.component_type, vid.component_index)
    #     i, j = node_map[string(comp["node_fr"])], node_map[string(comp["node_to"])]

    #     if val >= 0.5 !LightGraphs.has_edge(graph_nw, i, j)
    #         # Add the component back into the graph.
    #         LightGraphs.add_edge!(graph_nw, i, j)
    #         comp_sym = Dict(Symbol(key) => value for (key, value) in comp)
    #         comp_sym[:component_type] = vid.component_type
    #         MetaGraphs.set_props!(graph_nw, i, j, comp_sym)
    #     elseif val < 0.5 && LightGraphs.has_edge(graph_nw, i, j)
    #         # Remove the component from the graph.
    #         LightGraphs.rem_edge!(graph_nw, i, j)
    #     end
    # end

    # for component_type in ["pipe", "short_pipe"]
    #     for comp in values(wm_data[component_type])
    #         i, j = node_map[string(comp["node_fr"])], node_map[string(comp["node_to"])]
    #         q_sym = Symbol("q_" * component_type)
    #         q_val = JuMP.value(WM.var(wm, q_sym, comp["index"]))

    #         if LightGraphs.has_edge(graph_nw, i, j) && q_val < -1.0e-6
    #             LightGraphs.rem_edge!(graph_nw, i, j) 
    #             LightGraphs.add_edge!(graph_nw, j, i)

    #             comp_sym = Dict(Symbol(key) => value for (key, value) in comp)
    #             comp_sym[:component_type] = Symbol(component_type)
    #             comp_sym[:flow_abs] = abs(q_val)
    #             MetaGraphs.set_props!(graph_nw, j, i, comp_sym)
    #         elseif LightGraphs.has_edge(graph_nw, j, i) && q_val > 1.0e-6
    #             LightGraphs.rem_edge!(graph_nw, j, i)
    #             LightGraphs.add_edge!(graph_nw, i, j)

    #             comp_sym = Dict(Symbol(key) => value for (key, value) in comp)
    #             comp_sym[:component_type] = Symbol(component_type)
    #             comp_sym[:flow_abs] = abs(q_val)
    #             MetaGraphs.set_props!(graph_nw, i, j, comp_sym)
    #         end
    #     end
    # end
end