#-----------------------------------------------------------------------------#
# State (structure of arrays, genes × cells) and the synchronous step
#
# Every layer at tick t+1 is computed only from layers at tick t, so the
# cascade chromatin → binding → RNA → protein → (regulation) has a one-tick
# delay per modality. Two buffers are swapped after each step.
#-----------------------------------------------------------------------------#

"""One full set of per-cell layers, each `n_genes × n_cells`."""
struct Layers
    A::Matrix{Float32}   # chromatin accessibility, 0..1
    S::Matrix{Bool}      # promoter ON/OFF (telegraph state)
    R::Matrix{Int16}     # mRNA molecules
    P::Matrix{Int16}     # protein molecules
end

Layers(G::Integer, N::Integer) =
    Layers(zeros(Float32, G, N), zeros(Bool, G, N), zeros(Int16, G, N), zeros(Int16, G, N))

function Base.copyto!(dst::Layers, src::Layers)
    copyto!(dst.A, src.A); copyto!(dst.S, src.S)
    copyto!(dst.R, src.R); copyto!(dst.P, src.P)
    return dst
end

const DEFAULT_CHUNK = 256

"""
    Simulation(net, kin, programs, program_of_cell; seed = 1, chunk = 256)

`program_of_cell[c]` is the index (into `programs`) of cell `c`'s program.
Cells are processed in fixed chunks of `chunk` cells, each with its own RNG
seeded from `(seed, chunk index)`, so results are identical for any number of
threads.
"""
mutable struct Simulation
    const net::RegNetwork
    const kin::Kinetics
    const programs::Programs
    const program::Vector{Int16}
    cur::Layers
    nxt::Layers
    const field::Matrix{Float32}   # n_ligands × n_cells, local ligand level
    const rngs::Vector{Xoshiro}
    const chunk::Int
    const seed::UInt64
    tick::Int
end

# Independent stream per (seed, chunk): integer seeds are hashed (SHA-256) by Random.
_chunk_rng(seed::UInt64, k::Integer) = Xoshiro((UInt128(seed) << 64) | UInt128(k))

function Simulation(net::RegNetwork, kin::Kinetics, programs::Programs,
                    program_of_cell::AbstractVector{<:Integer};
                    seed::Integer = 1, chunk::Integer = DEFAULT_CHUNK)
    G = n_genes(net)
    N = length(program_of_cell)
    length(kin.k_on) == G || throw(DimensionMismatch("kinetics are for $(length(kin.k_on)) genes, network has $G"))
    size(programs.A0, 1) == G || throw(DimensionMismatch("programs are for a different network"))
    all(p -> 1 <= p <= n_programs(programs), program_of_cell) ||
        throw(ArgumentError("program index out of range 1:$(n_programs(programs))"))
    chunk >= 1 || throw(ArgumentError("chunk must be ≥ 1"))
    prog = Int16.(program_of_cell)
    cur = Layers(G, N)
    cur.A .= programs.A0[:, prog]
    nxt = Layers(G, N)
    copyto!(nxt, cur)
    s = UInt64(seed)
    rngs = [_chunk_rng(s, k) for k in 1:cld(N, chunk)]
    return Simulation(net, kin, programs, prog, cur, nxt,
                      zeros(Float32, n_ligands(net), N), rngs, Int(chunk), s, 0)
end

"""Single-program convenience: every cell runs all modules."""
Simulation(net::RegNetwork, kin::Kinetics, n_cells::Integer = 1; kw...) =
    Simulation(net, kin, Programs(net), ones(Int16, n_cells); kw...)

n_cells(sim::Simulation) = length(sim.program)
n_genes(sim::Simulation) = n_genes(sim.net)

function Base.show(io::IO, sim::Simulation)
    print(io, "Simulation(", n_genes(sim), " regulated genes × ", n_cells(sim),
          " cells, ", n_programs(sim.programs), " programs, tick ", sim.tick, ")")
end

"""Reset all layers to the initial state and reseed the RNGs."""
function reset!(sim::Simulation; seed::Integer = sim.seed)
    fill!(sim.cur.S, false); fill!(sim.cur.R, 0); fill!(sim.cur.P, 0)
    sim.cur.A .= sim.programs.A0[:, sim.program]
    copyto!(sim.nxt, sim.cur)
    fill!(sim.field, 0f0)
    for k in eachindex(sim.rngs)
        sim.rngs[k] = _chunk_rng(UInt64(seed), k)
    end
    sim.tick = 0
    return sim
end
