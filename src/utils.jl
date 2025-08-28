using Distributions

function add_noise(
    center::Union{Int64,Float64},
    std::Float64)
    lower = center - 3*std
    upper = center + 3*std
    dist = truncated(Normal(center, std), lower, upper)
    return round(Int, rand(dist))
end




function density_decay(
    density::Float64,
    shift_set::Int64,
    damp::Float64)::Int
    if density < 1.0
        return 0
    end
    density = range(0, density)
    density = damp .^ density
    shift_set = shift_set .* density
    shift_set = sum(shift_set)
    return round(Int,shift_set)
end

function wave_decay(
    ecosystem::Vector,
    shift_set::Int64,
    wave_damp::Float64,
    density_damp::Float64)::Int
    wave = range(0, length(ecosystem) - 1)
    wave = wave_damp .^ wave
    density = Vector{Float64}(undef, length(ecosystem))
    for e in eachindex(ecosystem)
        density[e] = density_decay(ecosystem[e], shift_set, density_damp)
    end
    shift_set = density .* wave
    shift_set = sum(shift_set)
    return round(Int,shift_set)
end


function cyclic_permuations(regulators::Vector{Int64})::Matrix{Int64}
    reg_seq = collect(regulators)
    return hcat([circshift(reg_seq, k) for k in 0:2]...)
end


function check_field_value(s::Any, field::Symbol)
    if !isstructtype(typeof(s))
        throw(ArgumentError("Input must be a struct"))
    end
    if !hasfield(typeof(s), field)
        throw(ArgumentError("Field $field does not exist in struct"))
    end
    
    value = getfield(s, field)
    # throw error - user can just remove faulty symbol
    if isnothing(value)
        throw(ArgumentError("Field $field exists but is empty (nothing)"))
    end
    return value
end


function make_unique_keys(dict::Dict, keys::Vector{String})
    # Initialize output vector for new keys
    new_keys = String[]
    
    for key in keys
        # Check if key exists in dictionary
        if haskey(dict, key)
            # Try to match key with pattern like "domain_0.1" or "celltype_1.2"
            suffix_match = match(r"^(.*?_\d+)\.(\d+)$", key)
            if !isnothing(suffix_match)
                # Key has format like "domain_0.1"
                base, suffix = suffix_match.captures
                suffix_num = parse(Int, suffix) + 1
                new_key = "$(base).$(suffix_num)"
                while haskey(dict, new_key)
                    suffix_num += 1
                    new_key = "$(base).$(suffix_num)"
                end
            else
                # Try to match key with pattern like "domain_0" or "celltype_1"
                base_match = match(r"^(.*?)_(\d+)$", key)
                if !isnothing(base_match)
                    # Key has format like "domain_0"
                    base, num = base_match.captures
                    suffix = 1
                    new_key = "$(base)_$(num).$(suffix)"
                    while haskey(dict, new_key)
                        suffix += 1
                        new_key = "$(base)_$(num).$(suffix)"
                    end
                else
                    # Key doesn't match either pattern (e.g., "simplekey")
                    base = key
                    suffix = 1
                    new_key = "$(base).$(suffix)"
                    while haskey(dict, new_key)
                        suffix += 1
                        new_key = "$(base).$(suffix)"
                    end
                end
            end
            push!(new_keys, new_key)
        else
            # If key doesn't exist, use it as is
            push!(new_keys, key)
        end
    end
    
    return new_keys
end

function logistic_sampling(
    range::Tuple{Float64,Float64},
    s::Float64 = 10.0)
    # Parameters
    a, b = range
    μ = b
    s = (b - a) / s 
    dist = Truncated(Logistic(μ, s), a, b)
    sample = rand(dist)
    return sample
end
#-----------------------------------------------------------------------------#
# Struct manipulation functions
#-----------------------------------------------------------------------------#
function update_cell_info!(
    cells::Vector{CellState},
    what::Vector{T},
    field::Symbol)::Vector{CellState} where T
    @inbounds for i in eachindex(cells)
        setfield!(cells[i].cell_info, field, what[i])
    end
    return cells
end

# function update_grns!(cells::Vector{CellState},
#     grn_set::Dict{String,GRN})::Vector{CellState}
#     @inbounds for c in eachindex(cells)
#         cell_info = values(cells[c].cell_info)
#         grn_local = Dict(k=> grn_set[k] for k in cell_info if haskey(grn_set, k))
#     end
# end
