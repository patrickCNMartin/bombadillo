using Distributions
using DelaunayTriangulation
using Graphs                 
using LinearAlgebra         
using SparseArrays
using Debugger


mutable struct tissue
    coordinates::AbstractMatrix
    spatial_graph::SimpleGraph
    distances::Dict
    barcodes::AbstractVector
    cells::AbstractVector
end

"""
    generate_coordinates()
"""
function generate_coordinates(
    x_range::Tuple = (1,1000),
    y_range::Tuple = (1,1000),
    z_range = nothing,
    n_cells::Int = 5000)::AbstractMatrix
    coordinates = Matrix{Float64}(undef, n_cells, 3)
    #coordinates[:,1] .= string.("cell_", 1:n_cells)
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




function create_spatial_graph(coordinates)::Tuple{SimpleGraph,Dict,AbstractMatrix}
    # Create Delaunay triangulation
    tri = DelaunayTriangulation.triangulate(transpose(coordinates[:,1:2]), delete_ghosts=true)
    
    #Extract edges from the triangulation
    edges = Set{Tuple{Int, Int}}()
    for triangle in tri.triangles
        push!(edges, minmax(triangle[1], triangle[2]))
        push!(edges, minmax(triangle[2], triangle[3]))
        push!(edges, minmax(triangle[3], triangle[1]))
        
    end
    #@debug
    # Create a graph
    n = size(coordinates, 1)
    g = Graphs.SimpleGraph(n)
    
    # # Add edges to the graph with Euclidean distances as weights
    weights = Dict{Tuple{Int, Int}, Float64}()
    for (u, v) in edges
        add_edge!(g, u, v)
        # Calculate Euclidean distance for this edge
        dist = norm(coordinates[u] .- coordinates[v])
        weights[(u, v)] = dist
        weights[(v, u)] = dist
    end
    return g, weights, coordinates
end

function get_neighborhood_graph(
    graph::SimpleGraph,
    initial_index::Int;
    depth::Int = 3)::Set{Int}

    visited = Set{Int}([initial_index])
    current_level = Set{Int}([initial_index])

    for _ in 1:depth
        next_level = Set{Int}()
        for u in current_level
            for neighbor in neighbors(graph, u)
                if neighbor ∉ visited
                    push!(next_level, neighbor)
                    push!(visited, neighbor)
                end
            end
        end
        current_level = next_level
        isempty(current_level) && break
    end
    return visited
end


function add_circles(
    graph::SimpleGraph,
    weights::Dict,
    coordinates::AbstractMatrix;
    n_circles::Int = 5,
    add_layers::Int = 3,
    depth_range::Tuple = (1,4),
    as_hotspots::Bool =  false,
    add_name::String = "Cell_1")
    #-------------------------------------------------------------------------#
    # Get points that will serve as center for circles
    # If as_hotspots = true, then we only those will be return
    #-------------------------------------------------------------------------#
    init_cell = 1
    center_indices = rand(1:nrow(coordinates), n_circles)
    for ci in center_indices
        if !as_hotspots
            local_depth = rand(depth_range[1]:depth_range[2])
            neighborhood = get_neighborhood_graph(
                graph,
                ci,
                depth = local_depth)
            init_labels[collect(neighborhood)] = init_cell
            init_cell += 1
        end
    end
    
end