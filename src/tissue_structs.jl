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
    wave_diffusion::Union{Int64, Nothing} = nothing
    wave_damp::Union{Float64, Nothing} = nothing
    density_damp::Union{Float64, Nothing} = nothing
    noise::Union{Float64, Nothing} = nothing
end

function Base.show(io::IO, ::MIME"text/plain",x::BaseCellType)
    println(io,"BaseCellType Struct:")
end

Base.@kwdef mutable struct Tissue
    coordinates::Union{AbstractVector, Nothing} = nothing
    cell_mesh::Union{Set, Nothing} = nothing
    cell_graph::Union{Dict, Nothing} = nothing
    distances::Union{AbstractVector, Nothing} = nothing
    barcodes::Union{AbstractVector, Nothing} = nothing
    cell_labels::Union{Dict, Nothing} = nothing
    ecosystems::Union{AbstractVector, Nothing} = nothing
    base_cell_types::Union{BaseCellType, Nothing} = nothing
    mesh::Union{DelaunayTriangulation.Triangulation, Nothing} = nothing # TMP will remove later
end

function Base.show(io::IO, ::MIME"text/plain",x::Tissue)
    println(io,"Tissue Struct:")
end


