#function add_no_good_cut!(wm::GenericWaterModel, optimizer::JuMP.OptimizerFactory, q_sol::Dict{Int, Float64}, h_sol::Dict{Int, Float64}, n::Int=wm.cnw)
function add_no_good_cut!(wm::GenericWaterModel, optimizer::JuMP.OptimizerFactory, n::Int=wm.cnw)
    xr_ones = Array{JuMP.VariableRef, 1}()
    xr_zeros = Array{JuMP.VariableRef, 1}()

    for (a, link) in wm.ref[:nw][n][:links_ne]
        xr = JuMP.value.(wm.var[:nw][n][:x_res][a])
        x_res, r_id = findmax(JuMP.value.(wm.var[:nw][wm.cnw][:x_res][a]))
        zero_indices = setdiff(1:length(xr), [r_id])
        xr_ones = vcat(xr_ones, wm.var[:nw][n][:x_res][a][r_id])
        xr_zeros = vcat(xr_zeros, wm.var[:nw][n][:x_res][a][zero_indices])
    end

    no_good_rhs = length(wm.ref[:nw][n][:links_ne]) - 1
    JuMP.@constraint(wm.model, sum(xr_ones) - sum(xr_zeros) <= no_good_rhs)
end
