function constraint_volume_sum(wm::WM.AbstractWaterModel, n::Int, i::Int, tank_ids::Vector{Int})
    V_sum = WM.var(wm, n, :V_sum, i)
    V_sum_tanks = sum(WM.var.(Ref(wm), n, :V, tank_ids))
    c = JuMP.@constraint(wm.model, V_sum == V_sum_tanks)
    append!(WM.con(wm, n, :volume_sum)[i], [c])
end