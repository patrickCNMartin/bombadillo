module Bombadillo


# abstract type Tissue end
# abstract type BaseCellType end

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
include("visualization.jl")

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
    add_spheres
end