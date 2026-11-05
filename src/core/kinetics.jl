#-----------------------------------------------------------------------------#
# Per-gene kinetic parameters (one tick = one time unit, default ~1 min)
#
# Telegraph promoter (bursting):
#   OFF → ON  with probability  k_on · A · f_binding(inputs)
#   ON  → OFF with probability  k_off
# RNA:     R' = R − Binomial(R, δ_r) + Poisson(α·[ON] + leak·A)
# Protein: P' = P − Binomial(P, δ_p) + Poisson(β·R)
# Chromatin: A' = A + κ·(A0 · f_chromatin(inputs) − A)
#
# Mean burst size is α / k_off transcripts; burst frequency is ≈ k_on·A·f.
#-----------------------------------------------------------------------------#

struct Kinetics
    k_on::Vector{Float32}
    k_off::Vector{Float32}
    alpha::Vector{Float32}
    leak::Vector{Float32}
    delta_r::Vector{Float32}
    beta::Vector{Float32}
    delta_p::Vector{Float32}
    kappa::Vector{Float32}
end

_pergene(x::Real, G) = fill(Float32(x), G)
function _pergene(x::AbstractVector, G)
    length(x) == G || throw(DimensionMismatch("expected $G per-gene values, got $(length(x))"))
    return Float32.(x)
end

"""
    Kinetics(G; k_on=0.5, k_off=0.3, alpha=15, leak=0.02, delta_r=0.3,
             beta=1.0, delta_p=0.07, kappa=0.05)

Scalars are broadcast to all `G` genes; vectors give per-gene values. Defaults
are Elowitz-repressilator-like at one tick per minute: mRNA half-life ≈ 2 min,
protein half-life ≈ 10 min, mean burst size α/k_off = 50 transcripts. With
these, `repressilator_network()` (K = 40, n = 2) oscillates with a period of
≈ 120–140 ticks.
"""
function Kinetics(G::Integer; k_on = 0.5, k_off = 0.3, alpha = 15, leak = 0.02,
                  delta_r = 0.3, beta = 1.0, delta_p = 0.07, kappa = 0.05)
    k = Kinetics(_pergene(k_on, G), _pergene(k_off, G), _pergene(alpha, G),
                 _pergene(leak, G), _pergene(delta_r, G), _pergene(beta, G),
                 _pergene(delta_p, G), _pergene(kappa, G))
    for f in fieldnames(Kinetics)
        v = getfield(k, f)
        all(x -> x >= 0, v) || throw(ArgumentError("$f must be non-negative"))
        f in (:k_on, :k_off, :delta_r, :delta_p, :kappa) && !all(<=(1), v) &&
            throw(ArgumentError("$f is a per-tick probability and must be ≤ 1"))
    end
    return k
end
Kinetics(net::RegNetwork; kw...) = Kinetics(n_genes(net); kw...)
