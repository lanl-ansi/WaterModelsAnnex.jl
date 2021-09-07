"Adds head gain constraints for pumps in `NC` formulations."
function WM.constraint_on_off_pump_head_gain(
    wm::AbstractNCZModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    q_min_forward::Float64,
)
    # Gather pump flow, head gain, and status variables.
    q = WM.var(wm, n, :q_pump, a)
    h_i = WM.var(wm, n, :h, node_fr)
    h_j = WM.var(wm, n, :h, node_to)
    z = WM.var(wm, n, :z_pump, a)

    # Add constraint equating head gain with respect to the pump curve.
    pump = WM.ref(wm, n, :pump, a)
    coeffs = WM._calc_head_curve_coefficients_quadratic(pump)
    q_expr = coeffs[1] * q^2 + coeffs[2] * q + coeffs[3]
    c_1 = JuMP.@NLconstraint(wm.model, z * q_expr <= z * (h_j - h_i))
    c_2 = JuMP.@NLconstraint(wm.model, z * q_expr >= z * (h_j - h_i))

    # Append the :on_off_pump_head_gain constraint array.
    append!(WM.con(wm, n, :on_off_pump_head_gain)[a], [c_1, c_2])
end


function WM.constraint_on_off_pump_head(
    wm::AbstractNCZModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
)
end


function WM.constraint_on_off_pump_power(
    wm::AbstractNCZModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    # Gather pump flow, scaled power, and status variables.
    q = WM.var(wm, n, :q_pump, a)
    P = WM.var(wm, n, :P_pump, a)
    z = WM.var(wm, n, :z_pump, a)

    # Add constraint equating power with respect to the power curve.
    power_qa = WM._calc_pump_power_quadratic_approximation(wm, n, a, z)
    c_1 = JuMP.@NLconstraint(wm.model, z * power_qa(q) <= z * P)
    c_2 = JuMP.@NLconstraint(wm.model, z * power_qa(q) >= z * P)

    # Append the :on_off_pump_power constraint array.
    append!(WM.con(wm, n, :on_off_pump_power)[a], [c_1, c_2])
end


function WM.constraint_on_off_pump_power_custom(wm::AbstractNCZModel, n::Int, a::Int, power_fixed::Float64, power_variable::Float64)
    # Gather pump flow, scaled power, and status variables.
    q = WM.var(wm, n, :q_pump, a)
    P = WM.var(wm, n, :P_pump, a)
    z = WM.var(wm, n, :z_pump, a)
    
    # Add constraint equating power with respect to the linear power curve.
    c = JuMP.@constraint(wm.model, z * (power_fixed + power_variable * q) == z * P)

    # Append the :on_off_pump_power constraint array.
    append!(WM.con(wm, n, :on_off_pump_power)[a], [c])
end
