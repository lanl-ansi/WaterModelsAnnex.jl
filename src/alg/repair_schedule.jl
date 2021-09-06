function turn_on_random_at_nw(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting}, nw::Int)
    nw_settings = control_settings[nw]
    off_ids = findall(i -> round(nw_settings.vals[i]) == 0.0, 1:length(nw_settings.vals))
    off_vars = nw_settings.variable_indices[off_ids]

    if length(off_ids) > 0
        random_var = Random.shuffle(off_vars)[1]
        
        if random_var.component_type == :pump
            pump_groups = WM.ref(wm, :pump_group)
            var_index = random_var.component_index
            group_id = findfirst(x -> var_index in x["pump_indices"], pump_groups)
            
            if group_id !== nothing
                pump_group = WM.ref(wm, :pump_group, group_id)
                pump_group_ids = sort(collect(pump_group["pump_indices"]))
                off_pump_vars = filter(x -> x.component_type == :pump, off_vars)
                off_pump_group_vars = filter(x -> x.component_index in pump_group_ids, off_pump_vars)
                var_index = sort(off_pump_group_vars, by = x -> x.component_index)[1]
                control_id = findfirst(x -> x == var_index, nw_settings.variable_indices)
                nw_settings.vals[control_id] = 1.0
                return true
            else
                control_id = findfirst(x -> x == random_var, nw_settings.variable_indices)
                nw_settings.vals[control_id] = 1.0
                return true
            end
        else
            control_id = findfirst(x -> x == random_var, nw_settings.variable_indices)
            nw_settings.vals[control_id] = 1.0
            return true
        end
    else
        return false
    end
end


function turn_off_random_at_nw(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting}, nw::Int)
    nw_settings = control_settings[nw]
    on_ids = findall(i -> round(nw_settings.vals[i]) == 1.0, 1:length(nw_settings.vals))
    on_vars = nw_settings.variable_indices[on_ids]

    if length(on_ids) > 0
        random_var = Random.shuffle(on_vars)[1]
        
        if random_var.component_type == :pump
            pump_groups = WM.ref(wm, :pump_group)
            var_index = random_var.component_index
            group_id = findfirst(x -> var_index in x["pump_indices"], pump_groups)
            
            if group_id !== nothing
                pump_group = WM.ref(wm, :pump_group, group_id)
                pump_group_ids = sort(collect(pump_group["pump_indices"]))
                on_pump_vars = filter(x -> x.component_type == :pump, on_vars)
                on_pump_group_vars = filter(x -> x.component_index in pump_group_ids, on_pump_vars)
                var_index = sort(on_pump_group_vars, by = x -> x.component_index)[end]
                control_id = findfirst(x -> x == var_index, nw_settings.variable_indices)
                nw_settings.vals[control_id] = 0.0
                return true
            else
                control_id = findfirst(x -> x == random_var, nw_settings.variable_indices)
                nw_settings.vals[control_id] = 0.0
                return true
            end
        else
            control_id = findfirst(x -> x == random_var, nw_settings.variable_indices)
            nw_settings.vals[control_id] = 0.0
            return true
        end
    else
        return false
    end
end


function filter_tank_infeasibilities(wm::WM.AbstractWaterModel, infeasibilities::Dict)
    tank_nodes = Set([x["node"] for (i, x) in WM.ref(wm, :tank)])
    return filter(x -> x.first.component_type == :node &&
        x.first.component_index in tank_nodes, infeasibilities)
end


function filter_out_reservoir_infeasibilities(wm::WM.AbstractWaterModel, infeasibilities::Dict)
    reservoir_nodes = [x["node"] for (i, x) in WM.ref(wm, :reservoir)]
    return filter(x -> !(x.first.component_type == :node &&
        (x.first.component_index in reservoir_nodes)), infeasibilities)
end


function filter_node_infeasibilities(wm::WM.AbstractWaterModel, infeasibilities::Dict)
    node_indices = Set([x["index"] for (i, x) in WM.ref(wm, :node)])
    return filter(x -> x.first.component_type == :node &&
        x.first.component_index in node_indices, infeasibilities)
end


function get_nearest_pumps(wm::WM.AbstractWaterModel, variable_index::WM._VariableIndex)
    graph, node_map = create_graph(wm.data)
    node_id = string(variable_index.component_index)
    node_to_pump_dist = Dict{Int, Int}()

    for pump in values(WM.ref(wm, :pump))
        from, to = node_map[node_id], node_map[string(pump["node_to"])]
        yen_state = LightGraphs.yen_k_shortest_paths(graph, from, to)
        node_to_pump_dist[pump["index"]] = yen_state.dists[1]
    end

    minimum_distance = minimum(values(node_to_pump_dist))
    min_pumps = filter(x -> x.second == minimum_distance, node_to_pump_dist)
    return sort(collect(keys(min_pumps)))
end


function get_nearest_valves(wm::WM.AbstractWaterModel, variable_index::WM._VariableIndex)
    graph, node_map = create_graph(wm.data)
    node_id = string(variable_index.component_index)
    node_to_valve_dist = Dict{Int, Int}()

    for (i, valve) in WM.ref(wm, :valve)
        from, to = node_map[node_id], node_map[string(valve["node_to"])]
        yen_state = LightGraphs.yen_k_shortest_paths(graph, from, to)
        node_to_valve_dist[valve["index"]] = yen_state.dists[1]
    end

    minimum_distance = minimum(values(node_to_valve_dist))
    min_valves = filter(x -> x.second == minimum_distance, node_to_valve_dist)
    return sort(collect(keys(min_valves)))
end


function get_pump_group(wm::WM.AbstractWaterModel, index::Int)
    pump_groups = WM.ref(wm, :pump_group)
    group_id = findfirst(x -> index in x["pump_indices"], pump_groups)
    return group_id !== nothing ? pump_groups[group_id] : nothing
end


function get_last_active_pump(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting}, index::Int, nw::Int)
    pump_group = get_pump_group(wm, index)

    if pump_group !== nothing
        variable_indices = control_settings[nw].variable_indices
        pump_group_indices = sort(collect(pump_group["pump_indices"]))
        pump_control_ids = findall(x -> x.component_type == :pump &&
            x.component_index in pump_group_indices, variable_indices)
        pump_control_vals = control_settings[nw].vals[pump_control_ids]
        last_active_id = findlast(x -> round(x) == 1.0, pump_control_vals)

        if last_active_id === nothing
            return index
        else
            return variable_indices[pump_control_ids[last_active_id]].component_index
        end
    else
        return index
    end
end


function get_first_inactive_pump(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting}, index::Int, nw::Int)
    pump_group = get_pump_group(wm, index)

    if pump_group !== nothing
        variable_indices = control_settings[nw].variable_indices
        pump_group_indices = sort(collect(pump_group["pump_indices"]))
        pump_control_ids = findall(x -> x.component_type == :pump &&
            x.component_index in pump_group_indices, variable_indices)
        pump_control_vals = control_settings[nw].vals[pump_control_ids]
        first_inactive_id = findfirst(x -> round(x) == 0.0, pump_control_vals)

        if first_inactive_id === nothing
            return index
        else
            return variable_indices[pump_control_ids[first_inactive_id]].component_index
        end
    else
        return index
    end
end


function turn_on_pump(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting}, i::Int, nw::Int)
    pump_id = get_first_inactive_pump(wm, control_settings, i, nw)
    index = findfirst(x -> x.component_type == :pump &&
        x.component_index == pump_id, control_settings[nw].variable_indices)

    if round(control_settings[nw].vals[index]) == 0.0
        control_settings[nw].vals[index] = 1.0
        return true
    else
        return false
    end
end


function turn_off_pump(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting}, i::Int, nw::Int)
    pump_id = get_last_active_pump(wm, control_settings, i, nw)
    index = findfirst(x -> x.component_type == :pump &&
        x.component_index == pump_id, control_settings[nw].variable_indices)

    if round(control_settings[nw].vals[index]) == 1.0
        control_settings[nw].vals[index] = 0.0
        return true
    else
        return false
    end
end


function calc_head_infeasibility_sum(infeasibilities::Dict)
    head_infs = filter(x -> x.first.component_type == :node, infeasibilities)

    if length(head_infs) > 0
        return sum(values(head_infs))
    else
        return 0.0
    end
end


function get_edges_from_path(path)
    edges = []

    for i in 1:length(path)-1
        push!(edges, LightGraphs.Edge(path[i], path[i+1]))
    end

    return edges
end


function get_edge_candidates_from_node(wm::WM.AbstractWaterModel, node_index::Int, graph)
    # Decompose the graph tuple into its constituents.
    graph_nw, node_map = graph[1], graph[2]
    
    # Determine the graph source node indices.
    tank_node_ids = [tank["node"] for tank in values(WM.ref(wm, :tank))]
    reservoir_node_ids = [res["node"] for res in values(WM.ref(wm, :reservoir))]
    source_node_ids = sort(vcat(tank_node_ids, reservoir_node_ids))
    source_vertex_ids = [node_map[string(x)] for x in source_node_ids]

    # Determine paths from the source nodes to the node with the infeasibility.
    node_vertex_id = node_map[string(node_index)]
    paths = vcat([collect(LightGraphs.all_simple_paths(graph_nw, s,
        node_vertex_id)) for s in source_vertex_ids]...)
    edges = collect(Set(vcat([get_edges_from_path(x) for x in paths]...)))

    # Select the subset of edges where controls can be changed.
    edge_types = [MetaGraphs.get_prop(graph_nw, x, :component_type) for x in edges]
    ids_to_change = findall(x -> x in ["pump", "valve", "regulator"], edge_types)
    return collect(Set(edges[ids_to_change]))
end


function get_edge_candidates_from_edge(wm::WM.AbstractWaterModel, index::Int, component_type::Symbol, graph)
    # Decompose the graph tuple into its constituents.
    graph_nw, node_map = graph[1], graph[2]
    
    # Determine the graph source node indices.
    tank_node_ids = [tank["node"] for tank in values(WM.ref(wm, :tank))]
    reservoir_node_ids = [res["node"] for res in values(WM.ref(wm, :reservoir))]
    source_node_ids = sort(vcat(tank_node_ids, reservoir_node_ids))
    source_vertex_ids = [node_map[string(x)] for x in source_node_ids]

    # Determine paths from the source nodes to the node with the infeasibility.
    i = node_map[string(WM.ref(wm, component_type, index, "node_fr"))]
    j = node_map[string(WM.ref(wm, component_type, index, "node_to"))]

    paths_to_i = vcat([collect(LightGraphs.all_simple_paths(graph_nw, s, i)) for s in source_vertex_ids]...)
    paths_to_j = vcat([collect(LightGraphs.all_simple_paths(graph_nw, s, j)) for s in source_vertex_ids]...)

    edges_to_i = collect(Set(vcat([get_edges_from_path(x) for x in paths_to_i]...)))
    edges_to_j = collect(Set(vcat([get_edges_from_path(x) for x in paths_to_j]...)))
    edges = collect(Set(vcat(edges_to_i, edges_to_j)))

    # Select the subset of edges where controls can be changed.
    edge_types = [MetaGraphs.get_prop(graph_nw, x, :component_type) for x in edges]
    ids_to_change = findall(x -> x in ["pump", "valve", "regulator"], edge_types)
    return collect(Set(edges[ids_to_change]))
end


function try_repair_node_max(wm::WM.AbstractWaterModel, node_index::Int, graph, control_settings::Vector{ControlSetting}, nw::Int)
    # Get the list of edges whose status could be changed.
    edge_candidates = get_edge_candidates_from_node(wm, node_index, graph)
    length(edge_candidates) == 0 && (return false)

    # Choose a controllable edge at random whose status could be changed.
    edge_to_change = Random.rand(edge_candidates)
    comp_to_change = MetaGraphs.props(graph[1], edge_to_change)

    if comp_to_change[:component_type] == "pump"
        # If the edge that was selected is a pump, try to turn it off.
        return turn_off_pump(wm, control_settings, comp_to_change[:index], nw)
    else
        # If the edge that was selected is not a pump, change the status.
        z_sym = Symbol("z_" * comp_to_change[:component_type])
        vid = WM._VariableIndex(nw, Symbol(comp_to_change[:component_type]), z_sym, comp_to_change[:index])
        vids = control_settings[nw].variable_indices
        control_id = findfirst(x -> compare_variable_indices(x, vid), vids)
        control_val = round(control_settings[nw].vals[control_id])
        control_settings[nw].vals[control_id] = control_val == 1.0 ? 0.0 : 1.0
        return true
    end
end


function try_repair_node_min(wm::WM.AbstractWaterModel, node_index::Int, graph, control_settings::Vector{ControlSetting}, nw::Int)
    # Get the list of edges whose status could be changed.
    edges_to_change = get_edge_candidates_from_node(wm, node_index, graph)
    length(edges_to_change) == 0 && (return false)

    # Choose a controllable edge at random whose status could be changed.
    edge_to_change = Random.rand(edges_to_change)
    comp_to_change = MetaGraphs.props(graph[1], edge_to_change)

    if comp_to_change[:component_type] == "pump"
        # If the edge that was selected is a pump, try to turn it on.
        return turn_on_pump(wm, control_settings, comp_to_change[:index], nw)
    else
        # If the edge that was selected is not a pump, change the status.
        z_sym = Symbol("z_" * comp_to_change[:component_type])
        vid = WM._VariableIndex(nw, Symbol(comp_to_change[:component_type]), z_sym, comp_to_change[:index])
        control_id = findfirst(x -> compare_variable_indices(x, vid), control_settings[nw].variable_indices)
        control_val = round(control_settings[nw].vals[control_id])
        control_settings[nw].vals[control_id] = control_val == 1.0 ? 0.0 : 1.0
        return true
    end
end


function try_repair_flow_max(wm::WM.AbstractWaterModel, index::Int, component_type::Symbol, graph, control_settings::Vector{ControlSetting}, nw::Int)
    # Select the subset of edges where controls can be changed.
    edges_to_change = get_edge_candidates_from_edge(wm, index, component_type, graph)
    length(edges_to_change) == 0 && (return false)

    # Choose an edge at random to change the control.
    edge_to_change = Random.shuffle(edges_to_change)[1]
    comp = MetaGraphs.props(graph[1], edge_to_change)
    comp_type, comp_index = comp[:component_type], comp[:index]

    # Get the index of the corresponding control variable.
    vid = WM._VariableIndex(nw, Symbol(comp_type), Symbol("z_" * comp_type), comp_index)
    variable_indices = control_settings[nw].variable_indices
    control_id = findfirst(x -> compare_variable_indices(x, vid), variable_indices)

    if round(control_settings[nw].vals[control_id]) == 1.0
        control_settings[nw].vals[control_id] = 0.0
        return true
    else
        control_settings[nw].vals[control_id] = 1.0
        return true
    end
end


function try_repair_flow_min(wm::WM.AbstractWaterModel, index::Int, component_type::Symbol, graph, control_settings::Vector{ControlSetting}, nw::Int)
    # Select the subset of edges where controls can be changed.
    edges_to_change = get_edge_candidates_from_edge(wm, index, component_type, graph)
    length(edges_to_change) == 0 && (return false)

    # Choose an edge at random to change the control.
    edge_to_change = Random.shuffle(edges_to_change)[1]
    comp = MetaGraphs.props(graph[1], edge_to_change)
    comp_type, comp_index = comp[:component_type], comp[:index]

    # Get the index of the corresponding control variable.
    vid = WM._VariableIndex(nw, Symbol(comp_type), Symbol("z_" * comp_type), comp_index)
    variable_indices = control_settings[nw].variable_indices
    control_id = findfirst(x -> compare_variable_indices(x, vid), variable_indices)

    if round(control_settings[nw].vals[control_id]) == 1.0
        control_settings[nw].vals[control_id] = 0.0
        return true
    else
        control_settings[nw].vals[control_id] = 1.0
        return true
    end
end


function repair_infeasibilities(wm::WM.AbstractWaterModel, graphs, control_settings::Vector{ControlSetting}, infeasibilities, nw_best::Int)
    nw_inf = sort([i.network_index for (i, x) in infeasibilities])[1]
    nw_ids = sort(collect(keys(control_settings)))
    nw_first_id = findfirst(x -> x == nw_inf, nw_ids)
    nws_change = nw_ids[max(1, nw_best - 1):nw_first_id]

    # Filter all tank infeasibilities to begin, if they exist.
    infeasibilities = filter_out_reservoir_infeasibilities(wm, infeasibilities)
    tank_infeasibilities = filter_tank_infeasibilities(wm, infeasibilities)
    node_infeasibilities = filter_node_infeasibilities(wm, infeasibilities)

    if length(tank_infeasibilities) > 0
        infeasibilities = tank_infeasibilities
    elseif length(node_infeasibilities) > 0
        infeasibilities = node_infeasibilities
    end

    # Select one of the infeasibilities encountered at random.
    infeasibility = Random.rand(collect(infeasibilities))
    # infeasibilities_arr = collect([(i, x) for (i, x) in infeasibilities])
    # sort!(infeasibilities_arr, by = x -> abs(x[2]));
    # infeasibility = infeasibilities_arr[end]

    println(JuMP.termination_status(wm.model), " ", infeasibility)

    tank_node_ids = [x["node"] for (i, x) in WM.ref(wm, :tank)]
    infeasibility_is_tank = infeasibility[1].component_type == :node &&
        infeasibility[1].component_index in tank_node_ids

    if (infeasibility_is_tank && infeasibility[2] < 0.0) && nw_inf == nw_ids[end-1]
        nw_change = nw_best = Random.rand(nw_ids[1:nw_first_id])
        changed = try_repair_node_max(wm, infeasibility[1].component_index,
            graphs[nw_change], control_settings, nw_change) 
    elseif (infeasibility_is_tank && infeasibility[2] > 0.0) && nw_inf == nw_ids[end-1]
        nw_change = nw_best = Random.rand(nw_ids[1:nw_first_id])
        changed = try_repair_node_min(wm, infeasibility[1].component_index,
            graphs[nw_change], control_settings, nw_change)
    elseif infeasibility[1].component_type == :node && infeasibility[2] < 0.0
        # If the head at a node is too large, try to change a component status.
        changed, nw_back = false, 0

        while !changed
            # Try to change a component status.
            (length(nws_change) == nw_back || nw_back > 2) && break
            nw_change = nws_change[end-nw_back]
            changed = try_repair_node_max(wm, infeasibility[1].component_index,
                graphs[nw_change], control_settings, nw_change)
            nw_back += 1
        end
    elseif infeasibility[1].component_type == :node && infeasibility[2] > 0.0
        # If the head at a node is too small, try to force a component on.
        changed, nw_back = false, 0

        while !changed
            # Try to change a component status.
            (length(nws_change) == nw_back || nw_back > 2) && break
            nw_change = nws_change[end-nw_back]
            changed = try_repair_node_min(wm, infeasibility[1].component_index,
                graphs[nw_change], control_settings, nw_change)
            nw_back += 1
        end
    elseif infeasibility[2] < 0.0
        # If the flow at an edge is too large, try to change a component status.
        changed, nw_back = false, 0

        while !changed
            # Try to change a component status.
            (length(nws_change) == nw_back || nw_back > 2) && break
            nw_change = nws_change[end-nw_back]
            changed = try_repair_flow_max(wm, infeasibility[1].component_index,
                infeasibility[1].component_type, graphs[nw_change], control_settings, nw_change)
            nw_back += 1
        end
    elseif infeasibility[2] > 0.0
        # If the flow at an edge is too small, try to change a component status.
        changed, nw_back = false, 0

        while !changed
            # Try to change a component status.
            (length(nws_change) == nw_back || nw_back > 2) && break
            nw_change = nws_change[end-nw_back]
            changed = try_repair_flow_min(wm, infeasibility[1].component_index,
                infeasibility[1].component_type, graphs[nw_change], control_settings, nw_change)
            nw_back += 1
        end
    end

    return nw_best
end


function repair_infeasibilities(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting}, infeasibilities, nw_best::Int)
    nw_inf = sort([i.network_index for (i, x) in infeasibilities])[1]
    nw_ids = sort(collect(keys(control_settings)))
    nw_first_id = findfirst(x -> x == nw_inf, nw_ids)
    nws_change = nw_ids[nw_best:max(nw_best, nw_first_id)]

    # Filter all tank infeasibilities to begin, if they exist.
    tank_infeasibilities = filter_tank_infeasibilities(wm, infeasibilities)
    infeasibilities = length(tank_infeasibilities) > 0 ? tank_infeasibilities : infeasibilities
    
    # Determine the time steps at which a control might be changed.
    nws_change = length(tank_infeasibilities) > 0 ? nw_ids[1:max(1, nw_first_id-1)] : nws_change

    # Select one of the infeasibilities encountered at random.
    infeasibility = Random.shuffle(collect(infeasibilities))[1]

    if infeasibility[1].component_type == :node && infeasibility[2] < 0.0
        # If the head at a node is too large, try to force a component off.
        turned_off, nw_back = false, 0

        while !turned_off
            # Try to turn off a random component.
            length(nws_change) == nw_back && break
            turned_off = turn_off_random_at_nw(wm,
                control_settings, nws_change[end-nw_back])
            nw_back += 1
        end
    elseif infeasibility[1].component_type == :node && infeasibility[2] > 0.0
        # If the head at a node is too small, try to force a component on.
        turned_on, nw_back = false, 0

        while !turned_on
            # Try to turn on a random component.
            turned_on = turn_on_random_at_nw(wm,
                control_settings, nws_change[end-nw_back])
            nw_back += 1
            length(nws_change) == nw_back && break
        end
    elseif infeasibility[2] > 0.0
        turn_off_random_at_nw(wm, control_settings, nws_change[end])
    elseif infeasibility[2] < 0.0
        turn_on_random_at_nw(wm, control_settings, nws_change[end])
    end
end


function repair_schedule(control_settings::Vector{ControlSetting}, network::Dict{String, <:Any}, nlp_optimizer)
    map(x -> x.vals = abs.(round.(x.vals)), control_settings)
    wm_cq = _instantiate_cq_model(network, nlp_optimizer)
    graphs = [create_graph(network) for nw in 1:network["time_series"]["num_steps"]]
    simulation_results, infeasibilities = simulate_control_settings_graph!(wm_cq, control_settings, graphs)
    nw_ids = sort(collect(keys(control_settings)))
    nw_best_id = nw_ids[1]

    iterations = 1

    while !all(x -> x.feasible, simulation_results) && iterations <= 50
        nw_inf = sort([i.network_index for (i, x) in infeasibilities])[1]
        nw_first_id = findfirst(x -> x == nw_inf, nw_ids)
        nw_best_id = max(nw_best_id, nw_first_id)
        
        nw_best_id = repair_infeasibilities(wm_cq, graphs, control_settings, infeasibilities, nw_best_id)
        simulation_results, infeasibilities = simulate_control_settings_graph!(wm_cq, control_settings, graphs)
        iterations += 1
    end

    feasible = all(x -> x.feasible, simulation_results)
    cost = sum(x.cost for x in simulation_results)
    println("Heuristic feasibility: $(feasible) --> Cost: $(cost)")
end


function repair_an_infeasibility(wm::WM.AbstractWaterModel, graphs, control_settings::Vector{ControlSetting}, infeasibilities)
    nw_inf = sort([i.network_index for (i, x) in infeasibilities])[1]
    nw_ids = sort(collect(keys(control_settings)))
    nw_first_id = findfirst(x -> x == nw_inf, nw_ids)
    nws_change = nw_ids[1:max(1, nw_first_id)]

    # Select one of the infeasibilities encountered at random.
    infeasibility = Random.rand(collect(infeasibilities))
    repair_infeasibility(wm, control_settings, infeasibility)
end


function repair_infeasibility(wm::WM.AbstractWaterModel, control_settings::Vector{ControlSetting}, infeasibility)
    nw_inf = infeasibility.first.network_index
    nw_ids = sort(collect(keys(control_settings)))
    nw_first_id = findfirst(x -> x == nw_inf, nw_ids)
    control_settings_new = deepcopy(control_settings)

    for nw in reverse(nw_ids[1:max(1, nw_first_id)])
        setting = control_settings_new[nw]

        for k in 1:length(setting.vals)
            setting.vals[k] = setting.vals[k] == 0.0 ? 1.0 : 0.0
            sim_res, infs = simulate_control_settings_sequential_detailed(wm, control_settings_new)

            if any(x -> compare_variable_indices(x.first, infeasibility.first), infs)
                setting.vals[k] = setting.vals[k] == 0.0 ? 1.0 : 0.0
            else
                return control_settings_new
            end
        end
    end

    return control_settings_new
end



# function repair_schedule(control_settings::Vector{ControlSetting}, network::Dict{String, <:Any}, nlp_optimizer)
#     map(x -> x.vals = abs.(round.(x.vals)), control_settings)
#     wm_cq = _instantiate_cq_model(network, nlp_optimizer)
#     simulation_results, infeasibilities = simulate_control_settings_sequential_detailed(wm_cq, control_settings)


#     # iterations = 1

#     # while !all(x -> x.feasible, simulation_results) && iterations <= 50
#     #     repair_infeasibilities(wm_cq, control_settings, infeasibilities)
#     #     simulation_results, infeasibilities = simulate_control_settings_sequential_detailed(wm_cq, control_settings)
#     #     infeasibilities === nothing && break
#     #     iterations += 1
#     # end

#     # feasible = all(x -> x.feasible, simulation_results)
#     # cost = sum(x.cost for x in simulation_results)
#     # println("Heuristic feasibility: $(feasible) --> Cost: $(cost)")
# end