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


#-----------------------------------------------------------------------------#
# New tissue generation
# Note: julia can handle fractions - could be nice to use them here
# simple tissue generator for now - we can add complexity to it later
#-----------------------------------------------------------------------------#
function generate_tissue(
    n_cells::Int64 = 5000,
    n_domains::Int64 = 5,
    n_types::Int64 = 10,
    n_genes::Int64 = 2000,
    coordinate_range::Tuple{Float64,Float64} = (0.0,1.0),
    max_diffusion::Float64 = 0.5,
    density_damp::Float64 = 0.1,
    diffusion_damp::Float64 = 0.3,
    static::Bool = true;
    bio_ref::Union{Nothing, BioRef} = nothing)
    #-------------------------------------------------------------------------#
    # For I will try a fully random just to make sure the code works
    #-------------------------------------------------------------------------#
    domain_vec = sample(1:n_domains, n_cells, replace = true)
    domain_vec = string.("domain_",string.(domain_vec))
    domain_labels = string.("domain_",string.(1:n_domains))
    type_vec = sample(1:n_types, n_cells, replace = true)
    type_vec = string.("celltype_",string.(type_vec))
    type_labels = string.("celltype_",string.(1:n_types))
    cell_vec = Vector{CellState}(undef,n_cells)
    #-------------------------------------------------------------------------#
    # Initialize Gene set - since there are a lot of other params
    # it might be useful to have the user set everything before hand
    #-------------------------------------------------------------------------#
    gene_set = initialize_genes(n_genes)
    #-------------------------------------------------------------------------#
    # Generate GRN sets for each cell type and domains
    #-------------------------------------------------------------------------#
    grn_set = Dict{String,GRN}()
    grns = vcat(domain_labels,type_labels)
    regulator_strength = gene_set.regulator_strength
    for g in eachindex(grns)
        grn = repressilator(regulator_strength)
        grn_set[grns[g]] = grn
    end
    #-------------------------------------------------------------------------#
    # Now we can generate a set of cells to populate the tissue
    # Note that we will collect the info to put into the TissueState struct
    # Just easier to store commonly used information rather than pulling 
    # it our each cell every time we might need it 
    #-------------------------------------------------------------------------#
    for c in 1:n_cells
        cell = initialize_cell(
            type_vec[c],
            domain_vec[c],
            grn_set,
            n_genes,
            coordinate_range)
        cell_vec[c] = cell
    end
    #-------------------------------------------------------------------------#
    # Pull and concat
    #-------------------------------------------------------------------------#
    coordinates = pull_coordinates(cell_vec)
    distances = cell_distances(coordinates, max_diffusion)
    #-------------------------------------------------------------------------#
    # build tissue
    #-------------------------------------------------------------------------#
    tissue = TissueState(
        cells = cell_vec,
        cell_types = type_vec,
        domains = domain_vec,
        coordinates = coordinates,
        cell_distances = distances,
        max_diffusion = max_diffusion,
        density_damp = density_damp,
        diffusion_damp = diffusion_damp,
        static = static)
    return tissue
end


function pull_coordinates(cell_vec::Vector{CellState})::Vector{Tuple{Float64,Float64,Float64}}
    coord = [c.coordinates for c in cell_vec]
    return coord
end



using SparseArrays
using NearestNeighbors
function cell_distances(coodinates::Vector{Tuple{Float64,Float64,Float64}},
    max_diffusion::Float64)::SparseMatrixCSC{Float64,Int64}
    #-------------------------------------------------------------------------#
    # Init and allocated
    #-------------------------------------------------------------------------#
    n = length(coodinates)
    tree = KDTree(hcat([[c[1], c[2], c[3]] for c in coodinates]...))
    I, J, V = Int[], Int[], Float64[]
    #-------------------------------------------------------------------------#
    # Loop over points - using this approach to avoid allocating 
    # unnessary points - don't care about cells that are too far away
    #-------------------------------------------------------------------------#
    for i in 1:n
        idxs, dists = inrange(tree,
            [coodinates[i][1], coodinates[i][2], coodinates[i][3]],
            max_diffusion)
        for (j, d) in zip(idxs, dists)
            if i < j
                push!(I, i, j)
                push!(J, j, i)
                push!(V, d, d)
            end
        end
    end
    sparse(I, J, V, n, n)
end
