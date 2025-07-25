using Distributions
using DelaunayTriangulation            
using LinearAlgebra         


function create_cell_mesh(
    base_tissue::Union{BaseTissue, Nothing} = nothing)::Tissue
    #-------------------------------------------------------------------------#
    # First create a tissue structure if not already done
    #-------------------------------------------------------------------------#
    if isnothing(base_tissue)
        base_tissue = BaseTissue()
    end
    tissue = Tissue()
    max_cells = base_tissue.max_cells
    min_angle = base_tissue.min_angle
    max_angle = base_tissue.max_angle
    max_area = base_tissue.max_area
    coordinate_range = base_tissue.coordinate_range
    init_factor = base_tissue.init_factor
    #-------------------------------------------------------------------------#
    # Generate coordinates using specified ranges
    # Create coordinates from intial random points 
    #-------------------------------------------------------------------------#
    coordinates = generate_coordinates_old(
        coordinate_range,
        max_cells,
        init_factor)
    tissue.seed_points = deepcopy(coordinates)
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
    niche, points, distances = spatial_graph(tri)

    #-------------------------------------------------------------------------#
    # rebuild
    #-------------------------------------------------------------------------#
    tissue.cell_graph = niche
    tissue.distances = distances
    tissue.coordinates = points
    tissue.barcodes = string.("cell_", 1:length(points))
    tissue.mesh = tri
    tissue.base_tissue = base_tissue
    return tissue
end

function generate_coordinates_old(
    coordinate_range::Tuple{Float64,Float64} = (0,1),
    max_cells::Int = 5000,
    init_factor::Float64 = 0.1)::AbstractVector
    n_cells = Int(ceil(max_cells * init_factor))
    x = rand(Uniform(coordinate_range[1],coordinate_range[2]), n_cells)
    y = rand(Uniform(coordinate_range[1],coordinate_range[2]), n_cells)
    coordinates = tuple.(x,y)
    return coordinates
end

function get_domains(base_tissue::BaseTissue)
    if base_tissue.domain_type isa Tuple{String,Int,Int}
        domains = [base_tissue.domain_type]
    else
        domains = base_tissue.domain_type
    end
    return domains
end


function add_domains(
    tissue::Tissue)::Tissue
    domains = get_domains(tissue.base_tissue)
    depth_range = tissue.base_tissue.depth_range
    init_labels = ones(length(tissue.cell_graph))
    graph = tissue.cell_graph
    for (d, n,l) in domains
        if d == "circle"
            init_labels = add_circles(
                init_labels,
                graph,
                n_circles = n,
                n_layers = l,
                depth_range = depth_range)
        else
            throw(DomainError("Unsupported domain type - Choose from => circle"))
        end
    end
    tissue.cell_labels = Int.(init_labels)
    return tissue
end


function add_circles(
    init_labels::Vector,
    graph::Vector;
    n_circles::Int = 5,
    n_layers::Int = 0,
    depth_range::Tuple{Int,Int} = (3,10),
    )::Vector 

    init_cell = maximum(init_labels) + 1
    center_indices = rand(1:length(init_labels), n_circles)
    for ci in center_indices
        local_depth = rand(depth_range[1]:depth_range[2])
        if n_layers == 0
            neighborhood, _ = get_neighborhood_graph(
                graph,
                ci,
                depth = local_depth)
            
            init_labels[neighborhood] .= init_cell
            init_cell += 1
        else
            if n_layers >= local_depth
                layer_depth = 1:local_depth
            else
                layer_depth = round.(Int, range(1, local_depth, length=n_layers+1))
            end
            prev_neighborhood = Set([ci])  
                
            for layer_idx in eachindex(layer_depth)
                current_depth = layer_depth[layer_idx]
                current_neighborhood, _ = get_neighborhood_graph(
                    graph,
                    ci,
                    depth = current_depth)
                    
                layer_nodes = setdiff(collect(current_neighborhood), collect(prev_neighborhood))
                
                
                init_labels[layer_nodes] .= init_cell
                prev_neighborhood = current_neighborhood
                init_cell += 1
            end
        end
    end
    return init_labels
end


