using Revise
using Debugger
using Pkg
Pkg.activate(".")
using Bombadillo
using BenchmarkTools
using MultivariateStats

# load test Data
sample = initialize_sample()
sample = add_spheres(sample)
sample = add_grns(sample,:domains)
sample = add_cells(sample,2,domain_bound = true)
sample = add_grns(sample,:cell_types)

checkpoint = deepcopy(sample)
checkpoint = add_counts(checkpoint)

gene_set = sample.grn_set["domain_5"].regulators
gene_set = string.("gene_",gene_set)

sample, pull = let_live(sample,50, gene_set, pull_layer=:rna_state)
# view_gene_cycle(sample, pull, merge_cells = true)

view_tissue_cycle(sample, [pull[1]],"test.gif",dim = 2)

# sample = add_counts(sample)

# export_sample(checkpoint, "tests/Data/checkpoint_spatial_sample")
# export_sample(sample, "tests/Data/cycle_spatial_sample")