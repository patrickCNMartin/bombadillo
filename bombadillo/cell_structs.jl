using DataFrame



"""
    cell

A mutable structure for each individual cell. We will check the cell type and
fetch the appropriate params from the cell_type struct to modulate the gene ranks
"""
mutable struct cell
    cell_id::UInt64
    type::AbstractString
    ecosystem::DataFrame
    gene_rank::AbstractVector
end

"""
    cell_type

For each cell type, we create a shift rule set which define where genes
shifts should occure for that cell type. 

This is essentially a param rule set. 
"""
struct cell_type
    type::AbstractString
    type_set::AbstractVector
    type_shift::AbstractVector
    contact_set::AbstractVector
    contact_shift::AbstractVector
    wave_set::AbstractVector
    wave_shift::AbstractVector
    wave_diffusion::Float64
    wave_damp::Float64
    noise::UInt16
end