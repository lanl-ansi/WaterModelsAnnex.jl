function constraint_flow_conservation_excess(
    wm::WM.AbstractWaterModel,
    i::Int;
    nw::Int = WM.nw_id_default,
)
    # Collect various indices for edge-type components connected to node `i`.
    pipe_fr, pipe_to = WM._get_from_and_to_components(wm, i, :pipe; nw = nw)
    des_pipe_fr, des_pipe_to = WM._get_from_and_to_components(wm, i, :des_pipe; nw = nw)
    pump_fr, pump_to = WM._get_from_and_to_components(wm, i, :pump; nw = nw)
    regulator_fr, regulator_to = WM._get_from_and_to_components(wm, i, :regulator; nw = nw)
    short_pipe_fr, short_pipe_to = WM._get_from_and_to_components(wm, i, :short_pipe; nw = nw)
    valve_fr, valve_to = WM._get_from_and_to_components(wm, i, :valve; nw = nw)

    # Collect various indices for node-type components connected to node `i`.
    reservoirs = WM.ref(wm, nw, :node_reservoir, i) # Reservoirs attached to node `i`.
    tanks = WM.ref(wm, nw, :node_tank, i) # Tanks attached to node `i`.
    demands = WM.ref(wm, nw, :node_demand, i) # Demands attached to node `i`.

    # Sum the constant demands required at node `i`.
    nondispatch_demand_ids = WM.ids(wm, nw, :nondispatchable_demand)
    nd_demands_at_i = filter(j -> j in nondispatch_demand_ids, demands)
    fixed_demands_at_i = WM.ref.(Ref(wm), nw, :nondispatchable_demand, nd_demands_at_i)
    fixed_flows_at_i = [x["flow_nominal"] for x in fixed_demands_at_i]    
    net_fixed_demand = length(fixed_flows_at_i) > 0 ? sum(fixed_flows_at_i) : 0.0

    # Get the indices of dispatchable demands connected to node `i`.
    dispatchable_demands = filter(j -> j in WM.ids(wm, nw, :dispatchable_demand), demands)

    # Initialize the flow conservation constraint dictionary entry.
    WM._initialize_con_dict(wm, :flow_conservation, nw = nw)

    # Add the flow conservation constraint.
    constraint_flow_conservation_excess(
        wm,
        nw,
        i,
        pipe_fr,
        pipe_to,
        des_pipe_fr,
        des_pipe_to,
        pump_fr,
        pump_to,
        regulator_fr,
        regulator_to,
        short_pipe_fr,
        short_pipe_to,
        valve_fr,
        valve_to,
        reservoirs,
        tanks,
        dispatchable_demands,
        net_fixed_demand,
    )
end


function constraint_pipe_flow_nonlinear(wm::WM.AbstractWaterModel, a::Int; nw::Int=WM.nw_id_default, kwargs...)
    node_fr, node_to = WM.ref(wm, nw, :pipe, a)["node_fr"], WM.ref(wm, nw, :pipe, a)["node_to"]
    exponent, L = WM.ref(wm, nw, :alpha), WM.ref(wm, nw, :pipe, a)["length"]
    base_length = get(wm.data, "base_length", 1.0)
    base_time = get(wm.data, "base_time", 1.0)

    r = WM._calc_pipe_resistance(WM.ref(wm, nw, :pipe, a), wm.data["head_loss"], wm.data["viscosity"], base_length, base_time)
    q_max_reverse = min(get(WM.ref(wm, nw, :pipe, a), "flow_max_reverse", 0.0), 0.0)
    q_min_forward = max(get(WM.ref(wm, nw, :pipe, a), "flow_min_forward", 0.0), 0.0)

    WM._initialize_con_dict(wm, :pipe_flow_nonlinear, nw=nw, is_array=true)
    WM.con(wm, nw, :pipe_flow_nonlinear)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_flow_nonlinear(wm, nw, a, node_fr, node_to, exponent, L, r, q_max_reverse, q_min_forward)
end


function constraint_volume_sum(wm::WM.AbstractWaterModel, i::Int; nw::Int=WM.nw_id_default, kwargs...)
    tank_ids = Vector{Int}(WM.ref(wm, nw, :tank_group, i, "tank_indices"))
    WM._initialize_con_dict(wm, :volume_sum, nw=nw, is_array=true)
    WM.con(wm, nw, :volume_sum)[i] = Array{JuMP.ConstraintRef}([])
    constraint_volume_sum(wm, nw, i, tank_ids)
end


function constraint_pipe_head_nonlinear(wm::WM.AbstractWaterModel, a::Int; nw::Int=WM.nw_id_default, kwargs...)
    node_fr, node_to = WM.ref(wm, nw, :pipe, a)["node_fr"], WM.ref(wm, nw, :pipe, a)["node_to"]
    exponent, L = WM.ref(wm, nw, :alpha), WM.ref(wm, nw, :pipe, a)["length"]
    base_length = get(wm.data, "base_length", 1.0)
    base_time = get(wm.data, "base_time", 1.0)

    r = WM._calc_pipe_resistance(WM.ref(wm, nw, :pipe, a), wm.data["head_loss"], wm.data["viscosity"], base_length, base_time)
    q_max_reverse = min(get(WM.ref(wm, nw, :pipe, a), "flow_max_reverse", 0.0), 0.0)
    q_min_forward = max(get(WM.ref(wm, nw, :pipe, a), "flow_min_forward", 0.0), 0.0)

    WM._initialize_con_dict(wm, :pipe_head_nonlinear, nw=nw, is_array=true)
    WM.con(wm, nw, :pipe_head_nonlinear)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_head_nonlinear(wm, nw, a, node_fr, node_to, exponent, L, r, q_max_reverse, q_min_forward)
end


function constraint_on_off_pump_flow_nonlinear(wm::WM.AbstractWaterModel, a::Int; nw::Int=WM.nw_id_default, kwargs...)
    node_fr, node_to = WM.ref(wm, nw, :pump, a)["node_fr"], WM.ref(wm, nw, :pump, a)["node_to"]
    q_min_forward = max(get(WM.ref(wm, nw, :pump, a), "flow_min_forward", 0.0), 0.0)
    coeffs =  WM._calc_head_curve_coefficients(WM.ref(wm, nw, :pump, a))

    WM._initialize_con_dict(wm, :on_off_pump_flow_nonlinear, nw=nw, is_array=true)
    WM.con(wm, nw, :on_off_pump_flow_nonlinear)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_flow_nonlinear(wm, nw, a, node_fr, node_to, coeffs, q_min_forward)
end


function constraint_on_off_pump_gain_nonlinear(wm::WM.AbstractWaterModel, a::Int; nw::Int=WM.nw_id_default, kwargs...)
    node_fr, node_to = WM.ref(wm, nw, :pump, a)["node_fr"], WM.ref(wm, nw, :pump, a)["node_to"]
    q_min_forward = max(get(WM.ref(wm, nw, :pump, a), "flow_min_forward", 0.0), 0.0)
    coeffs =  WM._calc_head_curve_coefficients(WM.ref(wm, nw, :pump, a))

    WM._initialize_con_dict(wm, :on_off_pump_gain_nonlinear, nw=nw, is_array=true)
    WM.con(wm, nw, :on_off_pump_gain_nonlinear)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_gain_nonlinear(wm, nw, a, node_fr, node_to, coeffs, q_min_forward)
end


function constraint_tank_nonlinear(wm::WM.AbstractWaterModel, i::Int; nw::Int=WM.nw_id_default, kwargs...)
    tank, node_index = WM.ref(wm, nw, :tank, i), WM.ref(wm, nw, :tank, i)["node"]
    WM._initialize_con_dict(wm, :tank_nonlinear, nw=nw, is_array=true)
    WM.con(wm, nw, :tank_nonlinear)[i] = Array{JuMP.ConstraintRef}([])
    constraint_tank_nonlinear(wm, nw, i, node_index)
end