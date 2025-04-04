using Distributions
using DelaunayTriangulation


"""
    generate_coordinates()
"""
function generate_coordinates(
    x_range = (1,1000),
    y_range = (1,1000),
    z_range = nothing,
    n_cells = 5000)::AbstractMatrix
    coordinates = Matrix{Float64}(undef, n_cells, 3)
    coordinates[:,1] .= rand(Uniform(x_range[1], x_range[2]), n_cells)
    coordinates[:,2] .= rand(Uniform(y_range[1], y_range[2]), n_cells)
    #-------------------------------------------------------------------------#
    # adding z dimension for now but not in use since we will have issues
    # with 3D triangulation 
    #-------------------------------------------------------------------------#
    if z_range !== nothing && z_range isa Tuple
        coordinates[:,3] .= rand(Uniform(z_range[1], z_range[2]), n_cells)
    else
        coordinates[:,3] .= ones(n_cells)
    end
    return coordinates
end


function spatial_graph(
    coordinates)
    triangles = triangulate(transpose(coordinates[:,1:2]))
    #spatial_nn = voronoi(triangles, clip = true, predicates = FastKernel())
    

end


function add_circles(
    coordinates,
    n_circles = 5,
    add_layers = 3,
    radial_range = (0.1,0.4))
    center_indices = rand(Uniform(1,nrow(coordinates)), n_circles)

end