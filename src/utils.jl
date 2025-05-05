using Distributions

function apply_noise(
    center::Int,
    std::Float64,
    n::Int)
    lower = center - 3*std
    upper = center + 3*std
    dist = truncated(Normal(center, std), lower, upper)
    return round(Int, rand(dist))
end


function count_unique_matrices(matvec::Vector{Matrix{T}}) where T
    counter = Dict{Tuple{Vararg{Tuple{Vararg{T}}}}, Int}()

    for mat in matvec
        key = Tuple.(eachrow(mat))  # each row → Tuple, whole matrix → Tuple of Tuples
        counter[key] = get(counter, key, 0) + 1
    end

    return counter
end


function damp_decay(initial_value::Float64, damp::Float64)
    vec = zeros(initial_value)
    factor = 1
    for i in 1:initial_value
        vec[i] = factor * damp 
        factor = factor * damp 
    end
    return vec
end