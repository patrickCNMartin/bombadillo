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


function add_domains_old(
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


#-----------------------------------------------------------------------------#
# New tissue generation
# Note: julia can handle fractions - could be nice to use them here
# simple tissue generator for now - we can add complexity to it later
#-----------------------------------------------------------------------------#
function initialize_tissue(
    n_cells::Int64 = 5000;
    tissue_dims::Int64 = 2,
    coordinate_range::Tuple{Float64,Float64} = (0.0,1.0),
    max_diffusion::Float64 = 0.5,
    density_damp::Float64 = 0.1,
    diffusion_damp::Float64 = 0.3,
    static::Bool = true)::TissueState
    #-------------------------------------------------------------------------#
    # Generate coordinates
    #-------------------------------------------------------------------------#
    x = generate_coordinates(coordinate_range, n_cells)
    y = generate_coordinates(coordinate_range, n_cells)
    if tissue_dims == 3
        z = generate_coordinates(coordinate_range,n_cells) 
    else
        z = zeros(length(x))
    end
    coordinates = tuple.(x,y,z)
    #-------------------------------------------------------------------------#
    # Cell distance matrix
    #-------------------------------------------------------------------------#
    distances = cell_distances(coordinates, max_diffusion)
    #-------------------------------------------------------------------------#
    # build tissue
    #-------------------------------------------------------------------------#
    tissue = TissueState(
        cell_types = nothing,
        domains = nothing,
        coordinates = coordinates,
        cell_distances = distances,
        max_diffusion = max_diffusion,
        density_damp = density_damp,
        diffusion_damp = diffusion_damp,
        static = static)
    
    return tissue
end



function generate_coordinates(
    coordinate_range::Tuple{Float64,Float64} = (0.0,1.0),
    n_cells::Int64 = 5000)
    return rand(Uniform(coordinate_range[1],coordinate_range[2]),n_cells)
end



function pull_coordinates(cell_vec::Vector{CellState})::Vector{Tuple{Float64,Float64,Float64}}
    coord = [c.coordinates for c in cell_vec]
    return coord
end


using SparseArrays
using NearestNeighbors
function cell_distances(coordinates::Vector{Tuple{Float64,Float64,Float64}},
    max_diffusion::Float64)::SparseMatrixCSC{Float64,Int64}
    #-------------------------------------------------------------------------#
    # Init and allocated
    #-------------------------------------------------------------------------#
    n = length(coordinates)
    mat = stack(coordinates)
    tree = KDTree(mat)
    I, J, V = Int[], Int[], Float64[]
    #-------------------------------------------------------------------------#
    # Loop over points - using this approach to avoid allocating 
    # unnessary points - don't care about cells that are too far away
    #-------------------------------------------------------------------------#
    @inbounds for i in 1:n
        ref_point = [coordinates[i][1], coordinates[i][2], coordinates[i][3]]
        idxs = inrange(tree, ref_point, max_diffusion, false)  # Indices only, unsorted
        for j in idxs
            if i < j  # Avoid duplicates and self-pairs
                dist = sqrt(sum((coordinates[i][k] - coordinates[j][k])^2 for k in 1:3))
                push!(I, i)
                push!(J, j)
                push!(V, dist)
                # push!(I, j)
                # push!(J, i)
                # push!(V, dist)  # Symmetric entry
            end
        end
    end
    sparse(I, J, V, n, n)
end
using NearestNeighbors
function reference_distances(coordinates::Vector{Tuple{Float64,Float64,Float64}},
    radius::Float64,
    reference_coordinates::Int64)
    #-------------------------------------------------------------------------#
    # Init and allocated
    #-------------------------------------------------------------------------#
    mat = stack(coordinates)
    tree = KDTree(mat)
    reference = stack(coordinates[reference_coordinates])
    idxs = inrange(tree,reference,radius, false)
    return idxs   
end

#-----------------------------------------------------------------------------#
# adding domains
# not sure where I want these to be called.
# In the generate_tissue function or should it be a seperate exported function
#-----------------------------------------------------------------------------#
# function add_domains(sample::SampleState,
#     n_domains::Vector{Int64} = [5],
#     domain_type::Vector{String} = ["sphere"])::SampleState
#     #-------------------------------------------------------------------------#
#     # Loop over domain types 
#     #-------------------------------------------------------------------------#
#     for d in domain_type

#     end
# end



function add_spheres(
    sample::SampleState,
    n_domains::Int64 = 5,
    domain_range::Tuple{Float64,Float64} = (0.1,0.3),
    allow_hollow = false)::SampleState
    #-------------------------------------------------------------------------#
    # Pull and check
    #-------------------------------------------------------------------------#
    min_range = domain_range[1]
    max_range = domain_range[2]
    if min_range <= 0.0
        throw(DomainError("Min sphere range should be > 0.0"))
    end
    coordinates = sample.tissue.coordinates
    domain_labels = sample.tissue.domains
    n_cells = sample.n_cells
    if isnothing(domain_labels)
        domain_labels = fill("domain_0", n_cells)
        low_bound = 1
    else 
        low_bound = [parse(Int, m.match) for m in match.(r"\d+$", domain_labels)]
        low_bound = maximum(low_bound) + 1
    end
    #-------------------------------------------------------------------------#
    # Loop - don't need to have a index here
    #-------------------------------------------------------------------------#
    for _ in 1:n_domains
        if !allow_hollow
            intial_coord = rand(1:length(coordinates))
            radius = rand(Uniform(min_range,max_range))
            sphere = reference_distances(coordinates,radius,intial_coord)
            domain = string("domain_",low_bound)
            domain_labels[sphere] .= domain
        else
            intial_coord = rand(1:length(coordinates))
            radius_1 = rand(Uniform(min_range,max_range))
            sphere_1 = reference_distances(coordinates,radius_1,intial_coord)
            radius_2 = rand(Uniform(min_range,max_range))
            sphere_2 = reference_distances(coordinates,radius_2,intial_coord)
            if radius_1 > radius_2
                sphere = setdiff(sphere_1,sphere_2)
            else
                sphere = setdiff(sphere_2,sphere_1)
            end
            domain = string("domain_",low_bound)
            domain_labels[sphere] .= domain
        end
        low_bound += 1
    end
    #-------------------------------------------------------------------------#
    # update sample, tissue, and cells 
    #-------------------------------------------------------------------------#
    sample.tissue.domains = domain_labels
    update_cell_info!(
        sample.cells,
        domain_labels,
        :domain)
    return sample
end


function add_spheres!(
    sample::SampleState,
    n_domains::Int64 = 5,
    domain_range::Tuple{Float64,Float64} = (0.1,0.3),
    allow_hollow = false)::SampleState
    #-------------------------------------------------------------------------#
    # Pull and check
    #-------------------------------------------------------------------------#
    min_range = domain_range[1]
    max_range = domain_range[2]
    if min_range <= 0.0
        throw(DomainError("Min sphere range should be > 0.0"))
    end
    coordinates = sample.tissue.coordinates
    domain_labels = sample.tissue.domains
    n_cells = sample.n_cells
    if isnothing(domain_labels)
        domain_labels = fill("domain_0", n_cells)
        low_bound = 1
    else 
        low_bound = [parse(Int, m.match) for m in match.(r"\d+$", domain_labels)]
        low_bound = maximum(low_bound) + 1
    end
    #-------------------------------------------------------------------------#
    # Loop - don't need to have a index here
    #-------------------------------------------------------------------------#
    for _ in 1:n_domains
        if !allow_hollow
            intial_coord = rand(1:length(coordinates))
            radius = rand(Uniform(min_range,max_range))
            sphere = reference_distances(coordinates,radius,intial_coord)
            domain = string("domain_",low_bound)
            domain_labels[sphere] .= domain
        else
            intial_coord = rand(1:length(coordinates))
            radius_1 = rand(Uniform(min_range,max_range))
            sphere_1 = reference_distances(coordinates,radius_1,intial_coord)
            radius_2 = rand(Uniform(min_range,max_range))
            sphere_2 = reference_distances(coordinates,radius_2,intial_coord)
            if radius_1 > radius_2
                sphere = setdiff(sphere_1,sphere_2)
            else
                sphere = setdiff(sphere_2,sphere_1)
            end
            domain = string("domain_",low_bound)
            domain_labels[sphere] .= domain
        end
        low_bound += 1
    end
    #-------------------------------------------------------------------------#
    # update sample, tissue, and cells 
    #-------------------------------------------------------------------------#
    sample.tissue.domains = domain_labels
    update_cell_info!(
        sample.cells,
        domain_labels,
        :domain)
    return sample
end