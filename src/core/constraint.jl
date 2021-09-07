
function constraint_flow_conservation_excess(
    wm::WM.AbstractWaterModel, n::Int, i::Int, pipe_fr::Array{Int64,1},
    pipe_to::Array{Int64,1}, des_pipe_fr::Array{Int64,1}, des_pipe_to::Array{Int64,1},
    pump_fr::Array{Int64,1}, pump_to::Array{Int64,1}, regulator_fr::Array{Int64,1},
    regulator_to::Array{Int64,1}, short_pipe_fr::Array{Int64,1},
    short_pipe_to::Array{Int64,1}, valve_fr::Array{Int64,1}, valve_to::Array{Int64,1},
    reservoirs::Array{Int64,1}, tanks::Array{Int64,1}, dispatchable_demands::Array{Int64,1},
    fixed_demand::Float64)
    # Collect flow variable references per component.
    q_pipe, q_des_pipe = WM.var(wm, n, :q_pipe), WM.var(wm, n, :q_des_pipe)
    q_pump, q_regulator = WM.var(wm, n, :q_pump), WM.var(wm, n, :q_regulator)
    q_short_pipe, q_valve = WM.var(wm, n, :q_short_pipe), WM.var(wm, n, :q_valve)
    q_reservoir, q_tank = WM.var(wm, n, :q_reservoir), WM.var(wm, n, :q_tank)
    q_demand, q_excess = WM.var(wm, n, :q_demand), WM.var(wm, n, :q_excess, i)

    # Add the flow conservation constraint.
    WM.con(wm, n, :flow_conservation)[i] = JuMP.@constraint(wm.model, -
        sum(q_pipe[a] for a in pipe_fr) + sum(q_pipe[a] for a in pipe_to) -
        sum(q_des_pipe[a] for a in des_pipe_fr) + sum(q_des_pipe[a] for a in des_pipe_to) -
        sum(q_pump[a] for a in pump_fr) + sum(q_pump[a] for a in pump_to) -
        sum(q_regulator[a] for a in regulator_fr) +
        sum(q_regulator[a] for a in regulator_to) -
        sum(q_short_pipe[a] for a in short_pipe_fr) +
        sum(q_short_pipe[a] for a in short_pipe_to) -
        sum(q_valve[a] for a in valve_fr) + sum(q_valve[a] for a in valve_to) == -
        sum(q_reservoir[k] for k in reservoirs) - sum(q_tank[k] for k in tanks) +
        sum(q_demand[k] for k in dispatchable_demands) + fixed_demand + q_excess)
end


function constraint_volume_sum(wm::WM.AbstractWaterModel, n::Int, i::Int, tank_ids::Vector{Int})
    V_sum = WM.var(wm, n, :V_sum, i)
    V_sum_tanks = sum(WM.var.(Ref(wm), n, :V, tank_ids))
    c = JuMP.@constraint(wm.model, V_sum == V_sum_tanks)
    append!(WM.con(wm, n, :volume_sum)[i], [c])
end