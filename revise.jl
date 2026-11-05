using Revise
using Debugger
using Pkg
Pkg.activate(".")
using Bombadillo
using BenchmarkTools
using MultivariateStats

using DataFrames

function abs_inplace!(dfs::Vector{DataFrame})
    for df in dfs
        for col in names(df)
            df[!, col] .= abs.(df[!, col])
        end
    end
    return dfs
end

# load test Data
sample = initialize_sample()
sample = add_spheres(sample)
sample = add_grns(sample,:domains)
# sample = add_cells(sample,2,domain_bound = true)
# sample = add_grns(sample,:cell_types)

checkpoint = deepcopy(sample)
checkpoint = add_counts(checkpoint)

gene_set = sample.grn_set["domain_5"].regulators
gene_set = string.("gene_",gene_set)

sample, pull = let_live(sample,1000, gene_set, pull_layer=:binding_state)
abs_inplace!(pull)
view_gene_cycle(sample, pull, merge_cells = true)

# view_tissue_cycle(sample, [pull[1]],"test.gif",dim = 2)

# sample = add_counts(sample)

# export_sample(checkpoint, "tests/Data/checkpoint_spatial_sample")
# export_sample(sample, "tests/Data/cycle_spatial_sample")