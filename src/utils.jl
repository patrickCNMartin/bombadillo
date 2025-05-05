using Distributions

function compute_shifts(
    center::Union{Int64,Float64},
    std::Float64)
    lower = center - 3*std
    upper = center + 3*std
    dist = truncated(Normal(center, std), lower, upper)
    return round(Int, rand(dist))
end


function ecosystem_cell_types(
    cells::Dict{Int, Tuple{Matrix{Float64}, Vector{Float64}}},
    cell_type::Dict)

    label_map = Dict(k => 1 for k in keys(cell_type))
    matrix_groups = Dict{UInt64, Vector{Int64}}() # hash of matrix => list of cell keys
    matrix_hash_map = Dict{UInt64, AbstractMatrix}()

    for (cell, (mat, _)) in cells
        h = hash(mat)
       
        if haskey(matrix_groups, h)
            if mat == matrix_hash_map[h]
                push!(matrix_groups[h], cell)
            else
                # Handle hash collision: rare but possible
                new_h = hash(rand(UInt))  # generate a new, unique hash
                while haskey(matrix_groups, new_h)
                    new_h = hash(rand(UInt))
                end
                matrix_groups[new_h] = [cell]
                matrix_hash_map[new_h] = mat
            end
        else
            matrix_groups[h] = [cell]
            matrix_hash_map[h] = mat
        end
    end
    

    # Assign labels
    for (i, cell_group) in enumerate(values(matrix_groups))
        for cell in cell_group
            label_map[cell] = i
        end
    end

    return label_map
end



function density_decay(initial_value::Float64, damp::Float64)
    initial_value = Int(initial_value)
    vec = zeros(initial_value)
    factor = 1
    for i in 1:initial_value
        vec[i] = factor * damp 
        factor = factor * damp 
    end
    return vec
end

function wave_decay(depth::Int64, damp::Float64)
    return damp^depth
end

