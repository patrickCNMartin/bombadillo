using CSV
using DataFrames

function export_tissue(
    tissue::Tissue,
    counts::Matrix{Float64},
    file_name = "spatial_data")
    coordinates = tissue.coordinates
    coordinates = DataFrame(
        barcodes = tissue.barcodes,
        x = first.(coordinates),
        y = last.(coordinates),
        cell_labels = tissue.cell_labels)
    coord_name = string(file_name,"_coordinates.csv")
    CSV.write(coord_name, coordinates)
    counts = DataFrame(counts,:auto)
    rename!(counts, Symbol.(tissue.barcodes))
    insertcols!(counts,1, :genes => string.("gene_",1:size(counts,1)))
    counts_name = string(file_name,"_counts.csv")
    CSV.write(counts_name,counts)
    return 0
end


function generate_sample()
    base_tissue = BaseTissue(domain_type = [("circle",5,0)])
    base_cells = BaseCells(
        type_shifts = (set_size = 50, max_shift = 1000, noise =  5.0),
        contact_shifts = (set_size = 50, max_shift = 1000, noise =  5.0),
        wave_shifts = (set_size = 50, max_shift = 1000, noise =  5.0),
        wave_damp = 0.5,
        density_damp = 0.01,
        wave_diffusion = 10)
    tissue = create_cell_mesh(base_tissue)
    tissue = add_domains(tissue)
    cells = add_cells(base_cells, tissue, shift_type = ["type_shifts","contact_shifts","wave_shifts"])
    tissue = add_ecosystems(tissue, cells)
    #tissue = add_ecotypes(tissue, collapse_to = 0)
    counts = add_genes(tissue, cells)
    export_tissue(tissue, counts, "../data/spatial_data_test")
    return tissue, cells, counts
end