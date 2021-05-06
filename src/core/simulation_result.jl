mutable struct SimulationResult
    feasible::Bool
    q_tank::Dict{Int, Float64}
    cost::Float64
end