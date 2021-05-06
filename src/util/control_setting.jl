mutable struct ControlSetting
    network_id::Int
    variable_indices::Array{WM._VariableIndex, 1}
    vals::Array{Float64, 1}
end