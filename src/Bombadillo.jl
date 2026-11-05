module Bombadillo

using Random: Random, AbstractRNG, Xoshiro

# New array-based core (genes × cells, Int16 counts, synchronous step)
include("core/sampling.jl")
include("core/network.jl")
include("core/kinetics.jl")
include("core/simulation.jl")
include("core/step.jl")
include("core/observe.jl")

export Edge, RegNetwork, Programs, Kinetics, Simulation, Trajectory,
    step!, run!, reset!, cell, set_level!, repressilator_network,
    n_genes, n_cells

# Legacy pipeline (to be ported onto the core, then removed)
include("tissue_structs.jl")
include("state_structs.jl")
include("graphs.jl")
include("utils.jl")
include("generate_cells.jl")
include("generate_tissue.jl")
include("generate_ecosystems.jl")
include("generate_genes.jl")
include("generate_grn.jl")
include("generate_samples.jl")
include("cycle.jl")

# Plotting lives in ext/BombadilloPlotsExt.jl and is loaded automatically
# when Plots is loaded alongside Bombadillo (`using Bombadillo, Plots`).
function plot_tissue end
function view_gene_cycle end
function view_tissue_cycle end

export Tissue,
    BaseTissue,
    BaseCells,
    create_cell_mesh,
    add_circles,
    add_cells,
    add_ecosystems,
    plot_tissue_2d,
    get_neighborhood_graph,
    compute_shifts,
    add_ecotypes,
    add_genes,
    export_tissue,
    generate_sample,
    generate_coordinates,
    add_domains,
    add_cells,
    add_ecotypes,
    repressilator,
    plot_tissue,
    initialize_sample,
    add_spheres,
    add_grns,
    initialize_cycling!,
    let_live,
    add_counts, 
    generate_gene_sample,
    export_sample,
    view_gene_cycle,
    view_tissue_cycle,
    delta_shift
end