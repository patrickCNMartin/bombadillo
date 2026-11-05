#-----------------------------------------------------------------------------#
# Running, recording and per-cell inspection
#-----------------------------------------------------------------------------#

const LAYERS = (:A, :S, :R, :P)

function _layer(sim::Simulation, layer::Symbol)
    layer in LAYERS || throw(ArgumentError("layer must be one of $(LAYERS), got :$layer"))
    return getfield(sim.cur, layer)
end

_gene_ids(sim::Simulation, ::Colon) = 1:n_genes(sim)
_gene_ids(sim::Simulation, g::Integer) = [Int(g)]
_gene_ids(sim::Simulation, g::AbstractString) = [sim.net.gene_index[g]]
_gene_ids(sim::Simulation, gs::AbstractVector) = reduce(vcat, (_gene_ids(sim, g) for g in gs); init = Int[])
_cell_ids(sim::Simulation, ::Colon) = 1:n_cells(sim)
_cell_ids(sim::Simulation, c::Integer) = [Int(c)]
_cell_ids(sim::Simulation, cs::AbstractVector{<:Integer}) = collect(Int, cs)

"""
    Trajectory

Recorded layer values, `data[gene, cell, frame]`, taken at `ticks[frame]`.
"""
struct Trajectory{T}
    layer::Symbol
    genes::Vector{String}
    cells::Vector{Int}
    ticks::Vector{Int}
    data::Array{T,3}
end

"""
    run!(sim, ticks; record = nothing, genes = :, cells = :, every = 1)

Advance `ticks` steps. With `record = :A | :S | :R | :P`, store that layer for
the selected genes and cells every `every` ticks (starting with the current
state) and return a [`Trajectory`](@ref); otherwise return `sim`.
"""
function run!(sim::Simulation, ticks::Integer; record::Union{Nothing,Symbol} = nothing,
              genes = Colon(), cells = Colon(), every::Integer = 1)
    every >= 1 || throw(ArgumentError("every must be ≥ 1"))
    if record === nothing
        for _ in 1:ticks
            step!(sim)
        end
        return sim
    end
    gi = _gene_ids(sim, genes)
    ci = _cell_ids(sim, cells)
    src = _layer(sim, record)
    nframes = ticks ÷ every + 1
    data = Array{eltype(src),3}(undef, length(gi), length(ci), nframes)
    tk = Vector{Int}(undef, nframes)
    _snapshot!(data, 1, _layer(sim, record), gi, ci)
    tk[1] = sim.tick
    f = 1
    for t in 1:ticks
        step!(sim)
        if t % every == 0
            f += 1
            _snapshot!(data, f, _layer(sim, record), gi, ci)
            tk[f] = sim.tick
        end
    end
    return Trajectory(record, sim.net.genes[gi], collect(Int, ci), tk, data)
end

function _snapshot!(data::Array{T,3}, f::Int, src::Matrix{T}, gi, ci) where {T}
    @inbounds for (j, c) in enumerate(ci), (i, g) in enumerate(gi)
        data[i, j, f] = src[g, c]
    end
    return data
end

"""
    cell(sim, c)

Inspect one cell: its program and every layer for its regulated genes, as a
NamedTuple of vectors (copies, safe to keep across steps).
"""
function cell(sim::Simulation, c::Integer)
    1 <= c <= n_cells(sim) || throw(BoundsError(sim.program, c))
    p = sim.program[c]
    return (; program = sim.programs.names[p], genes = sim.net.genes,
            A = sim.cur.A[:, c], S = sim.cur.S[:, c],
            R = sim.cur.R[:, c], P = sim.cur.P[:, c])
end

"""
    set_level!(sim, layer, gene, cells, value)

Paint a value into the current state, e.g. seed protein in a few cells:
`set_level!(sim, :P, "TetR", 1:10, 200)`.
"""
function set_level!(sim::Simulation, layer::Symbol, gene, cells, value)
    M = _layer(sim, layer)
    for g in _gene_ids(sim, gene), c in _cell_ids(sim, cells)
        M[g, c] = value
    end
    return sim
end

#-----------------------------------------------------------------------------#
# Reference circuits
#-----------------------------------------------------------------------------#

"""
    repressilator_network(; genes = ["lacI", "tetR", "cI"], K = 40, n = 2)

Three-gene ring of repressors (Elowitz & Leibler 2000): each gene represses
the next at the binding layer.
"""
function repressilator_network(; genes = ["lacI", "tetR", "cI"], K::Real = 40,
                               n::Real = 2, module_name = "repressilator")
    m = length(genes)
    m >= 3 && isodd(m) || throw(ArgumentError("a repressor ring needs an odd number ≥ 3 of genes"))
    edges = [Edge(genes[i], genes[mod1(i + 1, m)]; sign = -1, K, n, module_name)
             for i in 1:m]
    return RegNetwork(edges; genes)
end
