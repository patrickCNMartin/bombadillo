using DelaunayTriangulation
     

"""
    cell_type

For each cell type, we create a shift rule set which define where genes
shifts should occure for that cell type. 

This is essentially a param rule set. 
"""
Base.@kwdef mutable struct BaseCells
    types::Union{AbstractVector,Nothing} = nothing
    type_shifts::Union{Dict, NamedTuple} = (set_size = 10, max_shift = 500, noise = 10.0)
    contact_shifts::Union{Dict, NamedTuple} = (set_size = 10, max_shift = 500, noise = 10.0)
    wave_shifts::Union{Dict, NamedTuple} = (set_size = 10, max_shift = 500, noise = 10.0)
    chain_shifts::Union{Dict, NamedTuple} = (set_size = 10, max_shift = 500, noise = 10.0, chain=[1,2,3,4,5])
    shift_type::Vector{String} = ["type_shifts","contact_shifts","wave_shifts"]
    wave_diffusion::Int64 = 3
    wave_damp::Float64 = 0.1
    density_damp::Float64 = 0.1
    n_genes::Int64 = 2000
end

function Base.show(io::IO, ::MIME"text/plain",x::BaseCells)
    println(io,"BaseCells Struct:")
end

Base.@kwdef struct BaseTissue
    max_cells::Int64 = 5000
    coordinate_range::Tuple{Float64,Float64} = (0.0,1.0)
    init_factor::Float64 = 0.2
    min_angle::Float64 = 30.0
    max_angle::Float64 = 120.0
    max_area::Float64 = 0.01
    domain_type::Union{Tuple{String, Int64, Int64},Vector{Tuple{String, Int64, Int64}}} = ("circle", 5,0)
    depth_range::Tuple{Int64, Int64} = (5,10)

end

function Base.show(io::IO, ::MIME"text/plain",x::BaseTissue)
    println(io,"BaseTissue Struct:")
end



Base.@kwdef mutable struct Tissue
    base_tissue::Union{BaseTissue, Nothing} = nothing
    seed_points::Union{Vector, Nothing} = nothing
    coordinates::Union{Vector, Nothing} = nothing
    cell_graph::Union{Vector, Nothing} = nothing
    distances::Union{Vector, Nothing} = nothing
    barcodes::Union{Vector, Nothing} = nothing
    cell_labels::Union{Vector, Nothing} = nothing
    ecosystems::Union{Vector{Tuple{Matrix,Vector}}, Nothing} = nothing
    mesh::Union{DelaunayTriangulation.Triangulation, Nothing} = nothing # TMP will remove later
end

function Base.show(io::IO, ::MIME"text/plain",x::Tissue)
    println(io,"Tissue Struct:")
end

Base.@kwdef mutable struct Sample
    tissue::Union{Tissue, Nothing} = nothing
    counts::Union{Matrix, Nothing} = nothing
    base_cells::Union{BaseCells, Nothing} = nothing
    base_tissue::Union{BaseTissue, Nothing} = nothing
end
