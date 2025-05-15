

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
    #-------------------------------------------------------------------------#
    # relabel and repack
    #-------------------------------------------------------------------------#
    new_keys = zeros(maximum(vertices))
    new_keys[collect(vertices)] = 1:length(vertices)
    new_keys = Int.(new_keys)
    old_keys = collect(keys(niche))
    niches = Vector{Tuple{Int64, Set{Int64}}}(undef, length(old_keys))
    
    for idx in eachindex(niches)
        new_key = (idx, Set(new_keys[collect(niche[old_keys[idx]])]))
        niches[idx] = new_key
    end
    compact_dist = Vector{Tuple}(undef,length(edges))
    edges = collect(edges)
    for idx in eachindex(edges)
        e1 = new_keys[getindex(edges[idx],1)]
        e2 = new_keys[getindex(edges[idx],2)]
        compact_dist[idx] = (e1,e2,distances[idx])
    end
    return niches, points, compact_dist
end


function get_neighborhood_graph(
    graph::Vector{Tuple{Int64, Set{Int64}}},
    initial_index::Int;
    depth::Int = 3)::Tuple{Vector{Int64}, Vector{Int64}}

    # Internal Dict for fast lookup
    graph_dict = Dict(k => v for (k, v) in graph)

    current_level = Set([initial_index])
    visited_set = Set([initial_index])
    visited = [initial_index]
    distances = [0]

    for d in 1:depth
        next_level = Set{Int64}()
        for u in current_level
            neighbors = get(graph_dict, u, Set{Int64}())
            for neighbor in neighbors
                if neighbor ∉ visited_set
                    push!(visited_set, neighbor)
                    push!(visited, neighbor)
                    push!(distances, d)
                    push!(next_level, neighbor)
                end
            end
        end
        isempty(next_level) && break
        current_level = next_level
    end

    return visited, distances
end
