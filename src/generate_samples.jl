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
    base_tissue = BaseTissue(domain_type = [("circle",5,0)],
        init_factor=0.25)
    base_cells_type = BaseCells(
        type_shifts = (set_size = 50, max_shift = 1000, noise =  5.0),
        density_damp = 0.01)
    base_cells_contact = BaseCells(
        contact_shifts = (set_size = 50, max_shift = 1000, noise =  5.0),
        density_damp = 0.01)
    base_cells_wave = BaseCells(
        wave_shifts = (set_size = 50, max_shift = 1000, noise =  5.0),
        wave_damp = 0.5,
        density_damp = 0.01,
        wave_diffusion = 10)
    base_cells_all = BaseCells(
        type_shifts = (set_size = 50, max_shift = 1000, noise =  5.0),
        contact_shifts = (set_size = 50, max_shift = 1000, noise =  5.0),
        wave_shifts = (set_size = 50, max_shift = 1000, noise =  5.0),
        wave_damp = 0.5,
        density_damp = 0.01,
        wave_diffusion = 10)
    tissue = create_cell_mesh(base_tissue)
    tissue = add_domains(tissue)


    type = add_cells(base_cells_type, tissue, shift_type = ["type_shifts"])
    contact = add_cells(base_cells_contact,tissue,shift_type = ["contact_shifts"])
    wave = add_cells(base_cells_wave,tissue,shift_type = ["wave_shifts"])
    all = add_cells(base_cells_all,tissue,shift_type = ["type_shifts","contact_shifts","wave_shifts"])


    type_tissue = deepcopy(tissue)
    type_tissue = add_ecosystems(type_tissue, type)

    contact_tissue = deepcopy(tissue)
    contact_tissue = add_ecosystems(contact_tissue, contact)

    wave_tissue = deepcopy(tissue)
    wave_tissue = add_ecosystems(wave_tissue, wave)

    all_tissue = deepcopy(tissue)
    all_tissue = add_ecosystems(all_tissue, all)

    type_tissue = add_ecotypes(type_tissue, collapse_to = 0)
    contact_tissue = add_ecotypes(contact_tissue, collapse_to = 0)
    wave_tissue = add_ecotypes(wave_tissue, collapse_to = 0)
    all_tissue = add_ecotypes(all_tissue, collapse_to = 0)

    type_counts = add_genes(type_tissue, type)
    export_tissue(type_tissue, type_counts, "../data/spatial_data_type")

    contact_counts = add_genes(contact_tissue, contact)
    export_tissue(contact_tissue, contact_counts, "../data/spatial_data_contact")

    wave_counts = add_genes(wave_tissue, wave)
    export_tissue(wave_tissue, wave_counts, "../data/spatial_data_wave")

    all_counts = add_genes(all_tissue, all)
    export_tissue(all_tissue, all_counts, "../data/spatial_data_all")

    
    
    return 0
end

#-----------------------------------------------------------------------------#
# New sample generator
# Should add all the arguments for tissue, cell, and genes set here
#-----------------------------------------------------------------------------#
function initialize_sample(
    n_cells::Int64 = 5000, 
    n_genes::Int64 = 2000,
    batch::Int64 = 1;
    tissue_dims::Int64 = 2,
    coordinate_range::Tuple{Float64,Float64} = (0.0,1.0),
    max_diffusion::Float64 = 0.05,
    density_damp::Float64 = 0.1,
    diffusion_damp::Float64 = 0.3,
    static::Bool = true,
    max_leak::Int64 = 10,
    max_decay::Int64 = 100,
    translation_efficiency::Vector{Float64} = [0.5, 1.0])::SampleState
    #-------------------------------------------------------------------------#
    # Intialize tissue
    #-------------------------------------------------------------------------#
    tissue = initialize_tissue(
        n_cells,
        tissue_dims = tissue_dims,
        coordinate_range = coordinate_range,
        max_diffusion = max_diffusion,
        density_damp = density_damp,
        diffusion_damp = diffusion_damp,
        static = static)
    #-------------------------------------------------------------------------#
    # Initialize cells
    #-------------------------------------------------------------------------#
    cells = initialize_cells(
        tissue,
        n_genes,
        n_cells)
    #-------------------------------------------------------------------------#
    # Intialize Genes
    #-------------------------------------------------------------------------#
    genes = initialize_genes(
        n_genes,
        max_leak = max_leak,
        max_decay = max_decay,
        translation_efficiency = translation_efficiency)
    #-------------------------------------------------------------------------#
    # not much to do now but it will be useful later when we star moving to
    # atlas level sample generation
    #-------------------------------------------------------------------------#
    sample = SampleState(
        n_cells = n_cells,
        n_genes = n_genes,
        batch = batch,
        tissue = tissue,
        cells = cells,
        gene_state = genes)
    return sample
end