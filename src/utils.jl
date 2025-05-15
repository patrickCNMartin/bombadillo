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

