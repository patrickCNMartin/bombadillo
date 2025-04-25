
function add_ecosystems(
    tissue::Tissue;)::Tissue
    base_cells = tissue.base_cell_types
    wave_diffusion = base_cells.wave_diffusion
    wave_damp = base_cells.wave_damp
    ecosystems = Dict(k => Set{Int}() for k in keys(tissue.cell_graph))
    for eco in keys(ecosystems)
        local_eco = get_neighborhood_graph(
            tissue.cell_graph,
            eco,
            depth = wave_damp)
        local_eco = get_ecosystem(local_eco, tissue.cell_labels)
    end
end