

function spatial_graph(
    tri::Triangulation)
    #-------------------------------------------------------------------------#
    # Immutable struct so will need to pull out before filtering
    # can't do anything in place
    #-------------------------------------------------------------------------#
    edges = tri.graph.edges
    vertices = tri.graph.vertices
    niche = tri.graph.neighbours
    #-------------------------------------------------------------------------#
    # filter - remove negative indices
    # remove indices that are not in niche dict
    #-------------------------------------------------------------------------#
    for k in collect(keys(niche))
        k < 0 && delete!(niche, k)
    end
    keyset = Set(keys(niche))
    for s in values(niche)
        for x in copy(s)
            (x < 0 || !(x in keyset)) && delete!(s, x)
        end
    end

    filter!(t -> all(x -> x ≥ 0 && x in keys(niche), t), edges)
    filter!(s -> all(x -> x ≥ 0 && x in keys(niche), s), vertices)
    
    
    points = tri.points
    #-------------------------------------------------------------------------#
    # Compute distances
    #-------------------------------------------------------------------------#
    distances = Vector{Float64}(undef, length(edges))
    loc = 1
    for (u, v) in edges
        dist = norm(points[u] .- points[v])
        distances[loc] = dist
        loc += 1
    end
    points = [tri.points[i] for i in vertices]
    return edges, niche, points, distances
end

function get_neighborhood_graph(
    graph::Dict{Int64,Set{Int64}},
    initial_index::Int;
    depth::Int = 3)::Set{Int}

    visited = Set{Int}([initial_index])
    current_level = Set{Int}([initial_index])

    for _ in 1:depth
        next_level = Set{Int}()
        for u in current_level
            for neighbor in graph[u]
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
