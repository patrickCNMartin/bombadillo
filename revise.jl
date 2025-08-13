using Revise
using Debugger
using Pkg
Pkg.activate(".")
using Bombadillo
using BenchmarkTools


# load test Data
sample = initialize_sample()
sample = add_spheres(sample)
sample = add_grns(sample,:domains)
sample = add_cells(sample,2,domain_bound = true)
sample = add_grns(sample,:cell_types)