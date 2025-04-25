using Distributions
using DelaunayTriangulation            
using LinearAlgebra         


function create_cell_mesh(
    tissue::Union{Tissue, Nothing} = nothing,
    max_cells::Int = 5000,
    min_angle::Float64 = 30.0,
    max_angle::Float64 = 120.0,
    max_area::Float64 = 0.01,
    coordinate_range::Tuple{Int,Int} = (0,1),
    init_factor::Float64 = 0.2)::Tissue
    #-------------------------------------------------------------------------#
    # First create a tissue structure if not already done
    #-------------------------------------------------------------------------#
    if isnothing(tissue)
        tissue = Tissue()
    end
    #-------------------------------------------------------------------------#
    # Generate coordinates using specified ranges
    # Create coordinates from intial random points 
    #-------------------------------------------------------------------------#
    coordinates = generate_coordinates(
        coordinate_range,
        max_cells,
        init_factor)
    #-------------------------------------------------------------------------#
    # Now we generate triangles from the intial Set
    # and refine to the max number of cells and other params
    # Will likely generate data sets with different number of cells
    # Don't think that is a problem - I like the differential desinty
    #-------------------------------------------------------------------------#
    tri = DelaunayTriangulation.triangulate(coordinates;delete_ghosts = true)
    area = get_area(tri)
    refine!(tri;
        min_angle = min_angle,
        max_angle = max_angle,
        max_area = max_area * area,
        max_points = max_cells)
    #-------------------------------------------------------------------------#
    # Generate graph from mesh
    #-------------------------------------------------------------------------#
    edges, niche, points, distances = spatial_graph(tri)

    #-------------------------------------------------------------------------#
    # rebuild
    #-------------------------------------------------------------------------#
    tissue.cell_mesh = edges
    tissue.cell_graph = niche
    tissue.distances = distances
    tissue.coordinates = points
    tissue.barcodes = string.("cell_", 1:length(points))
    tissue.mesh = tri
    return tissue
end

function generate_coordinates(
    coordinate_range::Tuple{Int,Int} = (0,1),
    max_cells::Int = 5000,
    init_factor::Float64 = 0.1)::AbstractVector
    n_cells = Int(ceil(max_cells * init_factor))
    x = rand(Uniform(coordinate_range[1],coordinate_range[2]), n_cells)
    y = rand(Uniform(coordinate_range[1],coordinate_range[2]), n_cells)
    coordinates = tuple.(x,y)
    return coordinates
end




function add_circles(
    tissue::Tissue;
    n_circles::Int = 5,
    add_layers::Int = 0,
    depth_range::Tuple{Int,Int} = (3,10),
    as_hotspots::Bool =  false)::Tissue
    graph = tissue.cell_graph
    #-------------------------------------------------------------------------#
    # Get points that will serve as center for circles
    # If as_hotspots = true, then we only those will be return
    #-------------------------------------------------------------------------#
    init_labels = Dict(k => 1 for k in keys(graph))
    init_cell = 2
    center_indices = rand(keys(graph), n_circles)
    for ci in center_indices
        if !as_hotspots
            local_depth = rand(depth_range[1]:depth_range[2])
            if add_layers == 0
                neighborhood = get_neighborhood_graph(
                    graph,
                    ci,
                    depth = local_depth)
                foreach(k -> init_labels[k] = init_cell, neighborhood)
                init_cell += 1
            else
                if add_layers >= local_depth
                    layer_depth = 1:local_depth
                else
                    layer_depth = round.(Int, range(1, local_depth, length=add_layers+1))
                end
                prev_neighborhood = Set([ci])  
                
                for layer_idx in eachindex(layer_depth)
                    current_depth = layer_depth[layer_idx]
                    current_neighborhood = get_neighborhood_graph(
                        graph,
                        ci,
                        depth = current_depth)
                    
                    layer_nodes = setdiff(collect(current_neighborhood), collect(prev_neighborhood))
                
                    foreach(k -> init_labels[k] = init_cell, layer_nodes)
                    prev_neighborhood = current_neighborhood
                    init_cell += 1
                end
            end
        end
    end
    tissue.cell_labels = init_labels
    return tissue
end
