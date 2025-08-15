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

using CSV
using DataFrames

function export_sample(
    sample::SampleState,
    file_name = "spatial_data",
    labels::Symbol = :domains)
    coordinates = sample.tissue.coordinates
    barcodes = string.("cell_",collect(eachindex(coordinates)))
    coordinates = DataFrame(
        barcodes = barcodes,
        x = getindex.(coordinates,1),
        y = getindex.(coordinates,2),
        z = getindex.(coordinates,3),
        cell_labels = getfield(sample.tissue, labels))
    coord_name = string(file_name,"_coordinates.csv")
    CSV.write(coord_name, coordinates)
    counts = sample.out
    counts = DataFrame(counts,:auto)
    rename!(counts, Symbol.(barcodes))
    insertcols!(counts,1, :genes => string.("gene_",1:size(counts,1)))
    counts_name = string(file_name,"_counts.csv")
    CSV.write(counts_name,counts)
    return 0
end