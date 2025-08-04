using Revise
using Pkg
Pkg.activate(".")
using Bombadillo
using BenchmarkTools


# load test Data
sample = initialize_sample()
sample = add_spheres(sample)
sample = add_grns(sample,:domains)