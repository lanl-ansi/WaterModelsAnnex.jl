function degree_fr(data::Dict{String, <:Any}, component::Dict{String, <:Any})
    node_fr = component["node_fr"]

    pipes_fr = filter(x -> x.second["node_fr"] == node_fr, data["pipe"])
    pumps_fr = filter(x -> x.second["node_fr"] == node_fr, data["pump"])
    regulators_fr = filter(x -> x.second["node_fr"] == node_fr, data["regulator"])
    short_pipes_fr = filter(x -> x.second["node_fr"] == node_fr, data["short_pipe"])
    valves_fr = filter(x -> x.second["node_fr"] == node_fr, data["valve"])

    return length(pipes_fr) + length(short_pipes_fr)
    return length(pipes_fr) + length(pumps_fr) + length(regulators_fr) +
        length(short_pipes_fr) + length(valves_fr)
end


function degree_to(data::Dict{String, <:Any}, component::Dict{String, <:Any})
    node_to = component["node_to"]

    pipes_to = filter(x -> x.second["node_to"] == node_to, data["pipe"])
    pumps_to = filter(x -> x.second["node_to"] == node_to, data["pump"])
    regulators_to = filter(x -> x.second["node_to"] == node_to, data["regulator"])
    short_pipes_to = filter(x -> x.second["node_to"] == node_to, data["short_pipe"])
    valves_to = filter(x -> x.second["node_to"] == node_to, data["valve"])

    return length(pipes_to) + length(short_pipes_to)
    return length(pipes_to) + length(pumps_to) + length(regulators_to) +
        length(short_pipes_to) + length(valves_to)
end


function set_breakpoints_oa!(network_mn::Dict{String, <:Any}, result_mn::Dict{String, <:Any})
    for nw_id in sort([parse(Int, x) for x in collect(keys(network_mn["nw"]))])[1:end-1]
        network = network_mn["nw"][string(nw_id)]
        result = result_mn["solution"]["nw"][string(nw_id)]

        for (i, pipe) in network["pipe"]
            flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]

            if abs(result["pipe"][i]["q"]) > 1.0e-4
                flow_mid = max(flow_min, result["pipe"][i]["q"])
            else
                flow_mid = 0.5 * (flow_min + flow_max)
            end

            flow_min_mid = 0.5 * (flow_min + flow_mid)
            flow_max_mid = 0.5 * (flow_mid + flow_max)

            pipe["flow_lower_breakpoints"] = [flow_min,
                flow_min_mid, flow_mid, flow_max_mid, flow_max]
            pipe["flow_upper_breakpoints"] = [flow_min, flow_max]
         end

         for (i, pump) in network["pump"]
            flow_min, flow_max = pump["flow_min"], pump["flow_max"]
            
            if result["pump"][i]["q"] > 1.0e-4
                flow_mid = max(pump["flow_min"], result["pump"][i]["q"])
            else
                flow_mid = 0.5 * (flow_max + flow_min)
            end

            flow_min_mid = 0.5 * (flow_mid + flow_min)
            flow_max_mid = 0.5 * (flow_max + flow_mid)

            pump["flow_lower_breakpoints"] = [flow_min, flow_max]
            pump["flow_upper_breakpoints"] = range(flow_min, flow_max; length = 10)
        end
    end
end


function set_breakpoints_piecewise!(network_mn::Dict{String, <:Any}, result_mn::Dict{String, <:Any})
    for nw_id in sort([parse(Int, x) for x in collect(keys(network_mn["nw"]))])[1:end-1]
        network = network_mn["nw"][string(nw_id)]
        result = result_mn["solution"]["nw"][string(nw_id)]

        for (i, pipe) in network["pipe"]
            flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]

            if abs(result["pipe"][i]["q"]) > 1.0e-4
                flow_mid = max(flow_min, result["pipe"][i]["q"])
            else
                flow_mid = 0.5 * (flow_min + flow_max)
            end

            flow_min_mid = 0.5 * (flow_min + flow_mid)
            flow_max_mid = 0.5 * (flow_mid + flow_max)

            pipe["flow_lower_breakpoints"] = sort([flow_min,
                flow_min_mid, flow_mid, flow_max_mid, flow_max])
            pipe["flow_upper_breakpoints"] = sort([flow_min,
                flow_min_mid, flow_mid, flow_max_mid, flow_max])
         end

         for (i, pump) in network["pump"]
            flow_min, flow_max = pump["flow_min"], pump["flow_max"]            

            if result["pump"][i]["q"] > 1.0e-4
                flow_mid = max(pump["flow_min"], result["pump"][i]["q"])
            else
                flow_mid = 0.5 * (flow_max + flow_min)
            end

            flow_min_mid = 0.5 * (flow_mid + flow_min)
            flow_max_mid = 0.5 * (flow_max + flow_mid)

            pump["flow_lower_breakpoints"] = sort([flow_min,
                flow_min_mid, flow_mid, flow_max_mid, flow_max])
            pump["flow_upper_breakpoints"] = range(flow_min, flow_max; length = 10)
        end
    end
end


function set_breakpoints_piecewise_degree!(network_mn::Dict{String, <:Any}, result_mn::Dict{String, <:Any})
    for nw_id in sort([parse(Int, x) for x in collect(keys(network_mn["nw"]))])[1:end-1]
        network = network_mn["nw"][string(nw_id)]
        result = result_mn["solution"]["nw"][string(nw_id)]

        for (i, pipe) in network["pipe"]
            flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]

            if degree_fr(network, pipe) + degree_to(network, pipe) <= 2
                if abs(result["pipe"][i]["q"]) > 1.0e-4
                    flow_mid = max(flow_min, result["pipe"][i]["q"])
                else
                    flow_mid = 0.5 * (flow_min + flow_max)
                end

                flow_min_mid = 0.5 * (flow_mid + flow_min)
                flow_max_mid = 0.5 * (flow_max + flow_mid)

                pipe["flow_lower_breakpoints"] = [flow_min,
                    flow_min_mid, flow_mid, flow_max_mid, flow_max]
                pipe["flow_upper_breakpoints"] = [flow_min,
                    flow_min_mid, flow_mid, flow_max_mid, flow_max]
            else
                flow_mid = 0.5 * (flow_max + flow_min)
                flow_min_mid = 0.5 * (flow_mid + flow_min)
                flow_max_mid = 0.5 * (flow_max + flow_mid)

                pipe["flow_lower_breakpoints"] = [flow_min,
                    flow_min_mid, flow_mid, flow_max_mid, flow_max]
                pipe["flow_upper_breakpoints"] = [flow_min, flow_max]
            end
        end

        for (i, pump) in network["pump"]
            flow_min, flow_max = pump["flow_min_forward"], pump["flow_max"]

            if abs(result["pump"][i]["q"]) > flow_min
                flow_mid = result["pump"][i]["q"]
            else
                flow_mid = 0.5 * (flow_max + flow_min)
            end

            flow_min_mid = 0.5 * (flow_mid + flow_min)
            flow_max_mid = 0.5 * (flow_max + flow_mid)

            pump["flow_lower_breakpoints"] = [flow_min,
                flow_min_mid, flow_mid, flow_max_mid, flow_max]
            pump["flow_upper_breakpoints"] = range(flow_min, flow_max; length = 10)
        end
    end
end


function set_breakpoints_heuristic!(network_mn::Dict{String, <:Any}, result_mn::Dict{String, <:Any})
    for nw_id in sort([parse(Int, x) for x in collect(keys(network_mn["nw"]))])[1:end-1]
        network = network_mn["nw"][string(nw_id)]
        result = result_mn["solution"]["nw"][string(nw_id)]

        for (i, pipe) in network["pipe"]
            flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]
            pipe["flow_lower_breakpoints"] =
                range(flow_min, flow_max; length = 10)
            pipe["flow_upper_breakpoints"] =
                range(flow_min, flow_max; length = 10)

            # if degree_fr(network, pipe) + degree_to(network, pipe) <= 2
            #     if abs(result["pipe"][i]["q"]) > 1.0e-4
            #         flow_mid = max(flow_min, result["pipe"][i]["q"])
            #     else
            #         flow_mid = 0.5 * (flow_min + flow_max)
            #     end

            #     flow_min_mid = 0.5 * (flow_mid + flow_min)
            #     flow_max_mid = 0.5 * (flow_max + flow_mid)
                
            #     pipe["flow_lower_breakpoints"] =
            #         range(flow_min, flow_max; length = 25)
            #     pipe["flow_upper_breakpoints"] =
            #         range(flow_min, flow_max; length = 25)
            # else
            #     if abs(result["pipe"][i]["q"]) > 1.0e-4
            #         flow_mid = max(flow_min, result["pipe"][i]["q"])
            #     else
            #         flow_mid = 0.5 * (flow_min + flow_max)
            #     end

            #     flow_min_mid = 0.5 * (flow_mid + flow_min)
            #     flow_max_mid = 0.5 * (flow_max + flow_mid)

                
            #     pipe["flow_upper_breakpoints"] =
            #         range(flow_min, flow_max; length = 25)
            # end
        end

        for (i, pump) in network["pump"]
            flow_min, flow_max = pump["flow_min_forward"], pump["flow_max"]
            pump["flow_lower_breakpoints"] = 
                range(flow_min, flow_max; length = 10)
            pump["flow_upper_breakpoints"] =
                range(flow_min, flow_max; length = 10)

            # if abs(result["pump"][i]["q"]) > flow_min
            #     flow_mid = result["pump"][i]["q"]
            # else
            #     flow_mid = 0.5 * (flow_max + flow_min)
            # end

            # flow_min_mid = 0.5 * (flow_mid + flow_min)
            # flow_max_mid = 0.5 * (flow_max + flow_mid)            
        end
    end
end


function set_breakpoints_num!(network::Dict{String, <:Any}, num_breakpoints::Int)
    for pipe in values(network["pipe"])
        flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]
        breakpoints = range(flow_min, flow_max; length = num_breakpoints)
        pipe["flow_lower_breakpoints"] = breakpoints
        pipe["flow_upper_breakpoints"] = breakpoints
    end

    for pump in values(network["pump"])
        flow_min, flow_max = pump["flow_min_forward"], pump["flow_max"]
        breakpoints = range(flow_min, flow_max; length = num_breakpoints)
        pump["flow_lower_breakpoints"] = breakpoints
        pump["flow_upper_breakpoints"] = breakpoints
    end
end


function set_breakpoints_num_mn!(network::Dict{String, <:Any}, num_breakpoints::Int)
    for nw_id in sort([parse(Int, x) for x in collect(keys(network["nw"]))])[1:end-1]
        for pipe in values(network["nw"][string(nw_id)]["pipe"])
            flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]
            pipe["flow_lower_breakpoints"] = range(flow_min, flow_max; length = num_breakpoints)
            pipe["flow_upper_breakpoints"] = range(flow_min, flow_max; length = num_breakpoints)
        end
    
        for pump in values(network["nw"][string(nw_id)]["pump"])
            flow_min, flow_max = pump["flow_min_forward"], pump["flow_max"]
            pump["flow_lower_breakpoints"] = range(flow_min, flow_max; length = num_breakpoints)
            pump["flow_upper_breakpoints"] = range(flow_min, flow_max; length = num_breakpoints)
        end
    end
end