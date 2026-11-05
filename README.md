# Bombadillo

A Julia package for generating synthetic, procedurally built spatial omics data: a tissue
of cells, gene regulatory networks driving each cell's state, and counts sampled from that
state. Every dataset is reproducible from its parameters and seed.

Work in progress: an array-based simulation core (`src/core/`) is replacing the original
generators.

```julia
using Pkg; Pkg.develop(path = ".")
using Bombadillo
using Plots   # optional, enables the plotting functions
```

Tests: `julia --project=test test/runtests.jl`
