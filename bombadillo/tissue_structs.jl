using Distributions

"""
    generate_coordinates(max, min,)
"""
function generate_coordinates(
    x_range = (1,1000),
    y_range = (1,1000),
    z_range = nothing,
    n_cells = 5000)::AbstractMatrix
    coordinates = Matrix{Flaot64}(undef, n_cells, 3)
    coordinates[:,1] .= rand(Uniform(x_range[1], x_range[2]), n_cells)
    coordinates[:,2] .= rand(Uniform(y_range[1], y_range[2]), n_cells)
    if z_range !== nothing && z_range isa Tuple
        coordinates[:,3] .= rand(Uniform(z_range[1], z_range[2]), n_cells)
    else
        coordinates[:,3] .= ones(n_cells)
    end
    return coordinates
end