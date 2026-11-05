#-----------------------------------------------------------------------------#
# Regulatory network
#
# Genes 1:G are the *regulated* genes (anything that appears in an edge).
# Unregulated background genes are not simulated; they are generated at readout.
#
# Inputs to an edge are indexed in one space:
#   1:G        protein level of gene i in the same cell
#   G+1:G+L    local value of extracellular ligand field l (signalling)
#
# Edges act on one of two layers:
#   :chromatin  sets the target's accessibility A (slow, relaxes towards target)
#   :binding    sets the promoter's ON rate (fast, gated by A)
# A chromatin change at tick t therefore affects TF binding at t+1, RNA at t+2
# and protein at t+3: the modality cascade.
#
# Each edge belongs to a module (≤ 64). A cell's program is a set of modules,
# held as a UInt64 bitmask, so heterogeneous cells share one edge list.
#-----------------------------------------------------------------------------#

"""
    Edge(from, to; sign=-1, K=40, n=2, layer=:binding, module_name="core")

One regulatory interaction. `sign = +1` activates, `-1` represses. `K` is the
input level (molecules) at half-saturation and `n` the Hill coefficient.
`from` may be a gene or a ligand name.
"""
struct Edge
    from::String
    to::String
    sign::Int8
    K::Float32
    n::Float32
    layer::Symbol
    module_name::String
end

function Edge(from, to; sign::Integer = -1, K::Real = 40, n::Real = 2,
              layer::Symbol = :binding, module_name::AbstractString = "core")
    sign in (-1, 1) || throw(ArgumentError("sign must be +1 or -1, got $sign"))
    layer in (:binding, :chromatin) ||
        throw(ArgumentError("layer must be :binding or :chromatin, got :$layer"))
    K > 0 || throw(ArgumentError("K must be positive"))
    n > 0 || throw(ArgumentError("n must be positive"))
    return Edge(String(from), String(to), Int8(sign), Float32(K), Float32(n),
                layer, String(module_name))
end

"""
Edges of one layer, compressed by target gene (CSC): incoming edges of gene g
are `colptr[g]:colptr[g+1]-1`.
"""
struct EdgeSet
    colptr::Vector{Int32}
    src::Vector{Int32}
    K::Vector{Float32}
    n::Vector{Float32}
    sign::Vector{Int8}
    module_bit::Vector{UInt64}
end

n_edges(e::EdgeSet) = length(e.src)

struct RegNetwork
    genes::Vector{String}
    ligands::Vector{String}
    modules::Vector{String}
    chromatin::EdgeSet
    binding::EdgeSet
    gene_index::Dict{String,Int}
    module_index::Dict{String,Int}
    # module_genes[m] = genes touched (as source or target) by module m
    module_genes::Vector{Vector{Int}}
end

n_genes(net::RegNetwork) = length(net.genes)
n_ligands(net::RegNetwork) = length(net.ligands)
n_modules(net::RegNetwork) = length(net.modules)

"""
    RegNetwork(edges; genes = nothing, ligands = String[])

Build a network from a list of [`Edge`](@ref)s. Gene order is the order of
`genes` if given (genes with no edges are allowed: constitutive), otherwise
order of first appearance. Names in `ligands` are signalling inputs, not genes.
"""
function RegNetwork(edges::AbstractVector{Edge};
                    genes::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,
                    ligands::AbstractVector{<:AbstractString} = String[])
    ligs = String.(ligands)
    lig_set = Set(ligs)
    gnames = genes === nothing ? String[] : String.(genes)
    if genes === nothing
        seen = Set{String}()
        for e in edges, name in (e.from, e.to)
            (name in lig_set || name in seen) && continue
            push!(seen, name)
            push!(gnames, name)
        end
    end
    allunique(gnames) || throw(ArgumentError("duplicate gene names"))
    gidx = Dict(g => i for (i, g) in enumerate(gnames))
    G = length(gnames)

    mods = String[]
    for e in edges
        e.module_name in mods || push!(mods, e.module_name)
    end
    length(mods) <= 64 || throw(ArgumentError("at most 64 modules are supported"))
    midx = Dict(m => i for (i, m) in enumerate(mods))

    for e in edges
        e.to in lig_set && throw(ArgumentError("ligand $(e.to) cannot be a target"))
        haskey(gidx, e.to) || throw(ArgumentError("unknown target gene $(e.to)"))
        (haskey(gidx, e.from) || e.from in lig_set) ||
            throw(ArgumentError("unknown source $(e.from)"))
    end
    # one edge per (source, target, layer, module); duplicates are a model error
    keys_seen = Set{Tuple{String,String,Symbol,String}}()
    for e in edges
        k = (e.from, e.to, e.layer, e.module_name)
        k in keys_seen && throw(ArgumentError(
            "duplicate edge $(e.from) → $(e.to) ($(e.layer), $(e.module_name))"))
        push!(keys_seen, k)
    end

    src_index(name) = haskey(gidx, name) ? gidx[name] : G + findfirst(==(name), ligs)
    function edgeset(layer)
        sel = [e for e in edges if e.layer === layer]
        sort!(sel; by = e -> (gidx[e.to], src_index(e.from)))
        colptr = zeros(Int32, G + 1)
        colptr[1] = 1
        for e in sel
            colptr[gidx[e.to] + 1] += 1
        end
        cumsum!(colptr, colptr)
        # the kernel indexes with @inbounds, so the invariant is checked here
        colptr[end] == length(sel) + 1 || error("internal: malformed edge pointers")
        EdgeSet(colptr,
                Int32[src_index(e.from) for e in sel],
                Float32[e.K for e in sel],
                Float32[e.n for e in sel],
                Int8[e.sign for e in sel],
                UInt64[UInt64(1) << (midx[e.module_name] - 1) for e in sel])
    end

    module_genes = [Int[] for _ in mods]
    for e in edges
        mg = module_genes[midx[e.module_name]]
        for name in (e.from, e.to)
            haskey(gidx, name) && !(gidx[name] in mg) && push!(mg, gidx[name])
        end
    end
    foreach(sort!, module_genes)

    return RegNetwork(gnames, ligs, mods, edgeset(:chromatin), edgeset(:binding),
                      gidx, midx, module_genes)
end

function Base.show(io::IO, net::RegNetwork)
    print(io, "RegNetwork(", n_genes(net), " genes, ", n_ligands(net), " ligands, ",
          n_modules(net), " modules, ", n_edges(net.binding), " binding + ",
          n_edges(net.chromatin), " chromatin edges)")
end

#-----------------------------------------------------------------------------#
# Programs: which modules are active in a cell, and basal accessibility
#-----------------------------------------------------------------------------#

"""
    Programs(net, ["name" => ["module", ...], ...])

Maps a program name (a cell type, a domain, or a combination) to the modules
active in it. Genes touched by an active module start open (basal
accessibility 1); other regulated genes are closed (0) in that program. Edit
`A0` afterwards to paint partial accessibility.
"""
struct Programs
    names::Vector{String}
    mask::Vector{UInt64}
    A0::Matrix{Float32}    # G × n_programs basal accessibility
end

function Programs(net::RegNetwork, programs::AbstractVector{<:Pair})
    names = String[first(p) for p in programs]
    allunique(names) || throw(ArgumentError("duplicate program names"))
    mask = zeros(UInt64, length(programs))
    A0 = zeros(Float32, n_genes(net), length(programs))
    for (j, (_, mods)) in enumerate(programs)
        for m in mods
            haskey(net.module_index, m) || throw(ArgumentError("unknown module $m"))
            mi = net.module_index[m]
            mask[j] |= UInt64(1) << (mi - 1)
            A0[net.module_genes[mi], j] .= 1f0
        end
    end
    return Programs(names, mask, A0)
end
Programs(net::RegNetwork, programs::AbstractDict) = Programs(net, collect(programs))

"""All modules active, every gene open: the single-program case."""
Programs(net::RegNetwork) = Programs(net, ["all" => net.modules])

n_programs(p::Programs) = length(p.names)
