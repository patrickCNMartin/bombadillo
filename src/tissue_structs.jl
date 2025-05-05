using DelaunayTriangulation
     



"""
    cell_type

For each cell type, we create a shift rule set which define where genes
shifts should occure for that cell type. 

This is essentially a param rule set. 
"""
Base.@kwdef struct BaseCellType
    types::AbstractVector
    type_sets::Union{Dict, Nothing} = nothing
    type_shifts::Union{Dict, Nothing} = nothing
    contact_sets::Union{Dict, Nothing} = nothing
    contact_shifts::Union{Dict, Nothing} = nothing
    wave_sets::Union{Dict, Nothing} = nothing
    wave_shifts::Union{Dict, Nothing} = nothing
    wave_diffusion::Int64 = 5
    wave_damp::Float64 = 0.85
    density_damp::Float64 = 0.2
    noise::Float64 = 10.0
    n_genes::Int64 = 2000
end

function Base.show(io::IO, ::MIME"text/plain",x::BaseCellType)
    println(io,"BaseCellType Struct:")
end

Base.@kwdef mutable struct Cell
    type_shift::Union{AbstractVector, Nothing} = nothing
    contact_shift::Union{AbstractVector, Nothing} = nothing
    wave_shift::Union{AbstractVector, Nothing} = nothing
end

Base.@kwdef mutable struct Tissue
    coordinates::Union{AbstractVector, Nothing} = nothing
    cell_mesh::Union{Set, Nothing} = nothing
    cell_graph::Union{Dict, Nothing} = nothing
    distances::Union{AbstractVector, Nothing} = nothing
    barcodes::Union{AbstractVector, Nothing} = nothing
    cell_labels::Union{Dict, Nothing} = nothing
    ecosystems::Union{Dict, Nothing} = nothing
    base_cell_types::Union{BaseCellType, Nothing} = nothing
    mesh::Union{DelaunayTriangulation.Triangulation, Nothing} = nothing # TMP will remove later
end

function Base.show(io::IO, ::MIME"text/plain",x::Tissue)
    println(io,"Tissue Struct:")
end


