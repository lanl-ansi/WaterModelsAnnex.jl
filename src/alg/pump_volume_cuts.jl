function _collect_remaining_nws(wm::WM.AbstractWaterModel, n::Int)
    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(WM.nw_ids(wm)))
    n_id = findfirst(x -> x == n, network_ids)
    return network_ids[n_id:end]
end


function _sum_remaining_pump_flows(wm::WM.AbstractWaterModel, nws::Array{Int64, 1}, source_pumps::Array)
    return sum(sum(WM.var(wm, nw, :q_pump, a) for a in source_pumps) for nw in nws)
end


function _sum_remaining_reservoir_flows(wm::WM.AbstractWaterModel, nws::Array{Int64, 1})
    return sum(sum(WM.var(wm, nw, :q_reservoir, a) for a in WM.ids(wm, :reservoir)) for nw in nws)
end


function _sum_remaining_demands(wm::WM.AbstractWaterModel, nws::Array{Int64, 1})
    expr = WM.JuMP.AffExpr(0.0)

    for nw in nws
        expr += sum(WM.var(wm, nw, :q_demand)) # Sum over variable, dispatchable demands.
        expr += sum(x["flow_rate"] for (i, x) in WM.ref(wm, nw, :nondispatchable_demand))
    end

    return expr
end


function add_pump_volume_cuts!(wm::AbstractWaterModel)
    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(WM.nw_ids(wm)))

    # Start with the first network, representing the initial time step.
    n_1, n_f = network_ids[1], network_ids[end]
    source_pumps = [parse(Int, x) for x in find_source_pumps(wm.data["nw"]["1"])]

    for n in network_ids
        time_step = WM.ref(wm, n, :time_step)
        nws_remaining = _collect_remaining_nws(wm, n)
        demand_volume = _sum_remaining_demands(wm, nws_remaining)
        pump_volume = _sum_remaining_pump_flows(wm, nws_remaining, source_pumps)
        reservoir_volume = _sum_remaining_reservoir_flows(wm, nws_remaining)
        tank_volume = (sum(WM.var(wm, n_1, :V)) - sum(WM.var(wm, n, :V))) * inv(time_step)
        c = WM.JuMP.@constraint(wm.model, demand_volume + tank_volume <= pump_volume)
    end
end
