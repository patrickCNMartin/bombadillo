using Test
using Bombadillo
using Random, Statistics

const B = Bombadillo

# JET and AllocCheck compile a lot of code; skip them with BOMBADILLO_QUICK=1
const QUICK = get(ENV, "BOMBADILLO_QUICK", "0") == "1"

@testset "Bombadillo" begin
    include("core_sampling.jl")
    include("core_network.jl")
    include("core_dynamics.jl")
    QUICK || include("core_quality.jl")
end
