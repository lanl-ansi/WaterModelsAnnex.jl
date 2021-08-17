function add_tank_groups!(data::Dict{String,<:Any}, min_group_size::Int)
    data["tank_group"] = Dict{String, Any}()
    tank_ids = [x["index"] for (i, x) in data["tank"]]
    min_group_size = min(min_group_size, length(tank_ids))
    tank_group_id = 1 # The first group index    

    if length(tank_ids) > 1
        for k in min_group_size:length(tank_ids)
            for group in collect(Combinatorics.combinations(tank_ids, k))
                gid = string(tank_group_id)
                data["tank_group"][gid] = Dict{String, Any}("tank_indices" => group)
                tank_group_id += 1
            end
        end
    end
end