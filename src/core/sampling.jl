#-----------------------------------------------------------------------------#
# Allocation-free Poisson and binomial samplers for the hot loop.
#
# Distributions.jl picks a sampler type at runtime (a type-unstable branch),
# which allocates inside a per-cell, per-gene loop. These are the standard
# exact algorithms (Hörmann 1993) with a scalar API:
#   rand_poisson(rng, λ)    inversion for λ < 10, PTRS otherwise
#   rand_binomial(rng, n, p) inversion for n·min(p,1-p) < 10, BTRS otherwise
# Both return Int and never allocate.
#-----------------------------------------------------------------------------#

# log(k!) for k = 0:15 exactly, then the Stirling series (error < 1e-10).
const _LOGFACT_TABLE = ntuple(i -> sum(log, 1.0:Float64(i - 1); init = 0.0), 16)

@inline function logfactorial(k::Int)
    k < 16 && return @inbounds _LOGFACT_TABLE[k + 1]
    x = Float64(k + 1)
    ix2 = inv(x * x)
    return (x - 0.5) * log(x) - x + 0.9189385332046728 +
           (1.0 / 12 - (1.0 / 360 - (1.0 / 1260) * ix2) * ix2) / x
end

"""
    rand_poisson(rng, λ) -> Int

Exact Poisson draw. Inversion for `λ < 10`, Hörmann's PTRS otherwise.
"""
@inline function rand_poisson(rng::AbstractRNG, λ::Float64)::Int
    λ <= 0.0 && return 0
    if λ < 10.0
        # sequential inversion of the CDF
        p = exp(-λ)
        F = p
        u = rand(rng)
        k = 0
        while u > F && k < 200
            k += 1
            p *= λ / k
            F += p
        end
        return k
    end
    return _poisson_ptrs(rng, λ)
end
@inline rand_poisson(rng::AbstractRNG, λ::Real) = rand_poisson(rng, Float64(λ))

function _poisson_ptrs(rng::AbstractRNG, λ::Float64)::Int
    slam = sqrt(λ)
    loglam = log(λ)
    b = 0.931 + 2.53 * slam
    a = -0.059 + 0.02483 * b
    invalpha = 1.1239 + 1.1328 / (b - 3.4)
    vr = 0.9277 - 3.6224 / (b - 2.0)
    while true
        U = rand(rng) - 0.5
        V = rand(rng)
        us = 0.5 - abs(U)
        k = floor((2.0 * a / us + b) * U + λ + 0.43)
        if us >= 0.07 && V <= vr
            return Int(k)
        end
        (k < 0.0 || (us < 0.013 && V > us)) && continue
        ki = Int(k)
        if log(V) + log(invalpha) - log(a / (us * us) + b) <=
           -λ + ki * loglam - logfactorial(ki)
            return ki
        end
    end
end

"""
    rand_binomial(rng, n, p) -> Int

Exact binomial draw. Inversion when the expected count of the rarer outcome is
below 10, Hörmann's BTRS otherwise.
"""
@inline function rand_binomial(rng::AbstractRNG, n::Int, p::Float64)::Int
    (n <= 0 || p <= 0.0) && return 0
    p >= 1.0 && return n
    p > 0.5 && return n - rand_binomial(rng, n, 1.0 - p)
    n * p < 10.0 && return _binomial_inv(rng, n, p)
    return _binomial_btrs(rng, n, p)
end
@inline rand_binomial(rng::AbstractRNG, n::Integer, p::Real) = rand_binomial(rng, Int(n), Float64(p))

@inline function _binomial_inv(rng::AbstractRNG, n::Int, p::Float64)::Int
    q = 1.0 - p
    s = p / q
    a = (n + 1) * s
    r = q^n
    u = rand(rng)
    x = 0
    while u > r
        u -= r
        x += 1
        x > n && return n   # guards against round-off in the tail
        r *= (a / x - s)
    end
    return x
end

function _binomial_btrs(rng::AbstractRNG, n::Int, p::Float64)::Int
    q = 1.0 - p
    spq = sqrt(n * p * q)
    b = 1.15 + 2.53 * spq
    a = -0.0873 + 0.0248 * b + 0.01 * p
    c = n * p + 0.5
    vr = 0.92 - 4.2 / b
    alpha = (2.83 + 5.1 / b) * spq
    lpq = log(p / q)
    m = floor(Int, (n + 1) * p)
    h = logfactorial(m) + logfactorial(n - m)
    while true
        u = rand(rng) - 0.5
        v = rand(rng)
        us = 0.5 - abs(u)
        kf = floor((2.0 * a / us + b) * u + c)
        (kf < 0.0 || kf > n) && continue
        k = Int(kf)
        us >= 0.07 && v <= vr && return k
        v = log(v * alpha / (a / (us * us) + b))
        if v <= h - logfactorial(k) - logfactorial(n - k) + (k - m) * lpq
            return k
        end
    end
end
