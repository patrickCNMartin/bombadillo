module Bombadillo


# abstract type Tissue end
# abstract type BaseCellType end

include("tissue_structs.jl")
include("graphs.jl")
include("utils.jl")
include("generate_cells.jl")
include("generate_tissue.jl")
include("generate_ecosystems.jl")
include("visualization.jl")

export Tissue,
    BaseCellType,
    create_cell_mesh,
    add_circles,
    add_cells,
    add_ecosystems,
    plot_tissue_2d,
    get_neighborhood_graph,
    compute_shifts,
    ecosystem_cell_types
end