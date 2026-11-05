#-----------------------------------------------------------------------------#
# The kernel. `step_cells!` updates a contiguous range of cells with one RNG;
# it is allocation-free and type-stable (checked by the test suite).
#-----------------------------------------------------------------------------#

@inline function hill(x::Float32, K::Float32, n::Float32, s::Int8)
    r = (x / K)^n
    return s > Int8(0) ? r / (1f0 + r) : 1f0 / (1f0 + r)
end

"""
Combined effect of a gene's incoming edges in one layer: the product of the
Hill terms (AND logic) over edges whose module is active in the cell's program.
Returns 1 when no edge applies.
"""
@inline function regulation(es::EdgeSet, g::Int, mask::UInt64,
                            P::Matrix{Int16}, field::Matrix{Float32}, G::Int, c::Int)
    f = 1f0
    @inbounds for e in Int(es.colptr[g]):(Int(es.colptr[g + 1]) - 1)
        (es.module_bit[e] & mask) == 0 && continue
        s = Int(es.src[e])
        x = s <= G ? Float32(P[s, c]) : field[s - G, c]
        f *= hill(x, es.K[e], es.n[e], es.sign[e])
    end
    return f
end

@inline _sat16(x::Int) = (x > Int(typemax(Int16)) ? typemax(Int16) : x % Int16)

function step_cells!(nxt::Layers, cur::Layers, field::Matrix{Float32},
                     net::RegNetwork, kin::Kinetics, progs::Programs,
                     program::Vector{Int16}, rng::Xoshiro, cells::UnitRange{Int})
    G = size(cur.P, 1)
    chrom = net.chromatin
    bind = net.binding
    @inbounds for c in cells
        p = Int(program[c])
        mask = progs.mask[p]
        for g in 1:G
            A = cur.A[g, c]
            S = cur.S[g, c]
            R = Int(cur.R[g, c])
            P = Int(cur.P[g, c])

            # chromatin: relax towards the regulated target accessibility
            fc = regulation(chrom, g, mask, cur.P, field, G, c)
            nxt.A[g, c] = A + kin.kappa[g] * (progs.A0[g, p] * fc - A)

            # promoter: gated by last tick's accessibility
            fb = regulation(bind, g, mask, cur.P, field, G, c)
            u = rand(rng, Float32)
            nxt.S[g, c] = S ? (u >= kin.k_off[g]) : (u < kin.k_on[g] * A * fb)

            # RNA: from last tick's promoter state
            λr = (S ? kin.alpha[g] : 0f0) + kin.leak[g] * A
            Rn = R - rand_binomial(rng, R, Float64(kin.delta_r[g])) +
                 rand_poisson(rng, Float64(λr))
            # protein: from last tick's RNA
            Pn = P - rand_binomial(rng, P, Float64(kin.delta_p[g])) +
                 rand_poisson(rng, Float64(kin.beta[g]) * R)
            nxt.R[g, c] = _sat16(Rn)
            nxt.P[g, c] = _sat16(Pn)
        end
    end
    return nothing
end

_chunk_range(sim::Simulation, k::Int) =
    ((k - 1) * sim.chunk + 1):min(k * sim.chunk, n_cells(sim))

"""
    step!(sim)

Advance one tick. Chunks run on all available threads; results do not depend
on the thread count.
"""
function step!(sim::Simulation)
    nchunks = length(sim.rngs)
    if Threads.nthreads() == 1 || nchunks == 1
        for k in 1:nchunks
            step_cells!(sim.nxt, sim.cur, sim.field, sim.net, sim.kin, sim.programs,
                        sim.program, sim.rngs[k], _chunk_range(sim, k))
        end
    else
        Threads.@threads :dynamic for k in 1:nchunks
            step_cells!(sim.nxt, sim.cur, sim.field, sim.net, sim.kin, sim.programs,
                        sim.program, sim.rngs[k], _chunk_range(sim, k))
        end
    end
    sim.cur, sim.nxt = sim.nxt, sim.cur
    sim.tick += 1
    return sim
end
