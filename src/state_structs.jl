#-----------------------------------------------------------------------------#
# Define const
#-----------------------------------------------------------------------------#
const temporal_state = [1,2,3,4,5,6]
#-----------------------------------------------------------------------------#
# depends
#-----------------------------------------------------------------------------#
using SparseArrays
using Parameters
#-----------------------------------------------------------------------------#
# Gene struct - no need for mutability once set. 
#-----------------------------------------------------------------------------#
Base.@kwdef struct GeneState
    n_genes::Int64 = 2000
    genes::Union{Vector{Int64}, Vector{String}}
    saturation_rank::Vector{Int64}
    leak_rate::Vector{Int64} # if 0 then no change in baseline rank if more then it will go up till saturation of not repressed
    decay_rate::Vector{Int64}
    translation_efficiency::Vector{Float64}
    regulator_strength::Vector{Float64}
end
function Base.show(io::IO, ::MIME"text/plain",x::GeneState)
    println(io,"GeneState Struct")
end

#-----------------------------------------------------------------------------#
# GRN struct - we use this to define GRN types - complex regulatory patterns
# E.g how do we encode a reprisilator? Will need to modify this
# No need for mutability once set.
#-----------------------------------------------------------------------------#
Base.@kwdef struct GRN
    regulatory_rel::Matrix{Int64}# columns as GRN "steps"
    regulators::Vector{Int64} # What genes are involved?
    regulator_strength::Vector{Float64} # How strongly are they "on" - master regulator always 1
    messaging_output::Union{Nothing,Vector{Int64}} = nothing# Messages that will be diffused out
    metabolic_output::Union{Nothing,Vector{Int64}} = nothing# Catch all for everything else
    chromatin_remodelling::Union{Nothing,Vector{Int64}} = nothing # If they remodel chromatin - where do they do that?
    tf_binding::Union{Nothing,Vector{Int64}} = nothing# If they bind somewhere where do they do that?
end

function Base.show(io::IO, ::MIME"text/plain",x::GRN)
    println(io,"GRN Struct")
end
#-----------------------------------------------------------------------------#
# CellInfo - meta data storage for fast lookup and replacement
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct CellInfo
    celltype::Union{String, Nothing} = nothing
    domain::Union{String, Nothing} = nothing
end

#-----------------------------------------------------------------------------#
# Cell States struct
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct CellState
    cell_info::CellInfo = nothing # Maybe using a dict to keep all the meta info will be better
    cycle_position::Int64
    ecosystem::Union{Matrix{Float64},Nothing} = nothing
    coordinates::Tuple{Float64,Float64,Float64}
    chromatin_state::Union{SparseVector,Nothing} = nothing
    binding_state::Union{SparseVector,Nothing} = nothing
    rna_state::Union{Vector{Int64}, Nothing} = nothing
    protein_state::Union{Vector{Int64}, Nothing} = nothing
    metabolome_state::Union{SparseVector,Nothing} = nothing
    messaging_state::Union{SparseVector,Nothing} = nothing
end

function Base.show(io::IO, ::MIME"text/plain",x::CellState)
    println(io,"CellState Struct:")
end

#-----------------------------------------------------------------------------#
# Tissue struct - general tissue information 
# Starting state I guesss since we want to add dynamic cells and mvt
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct TissueState 
    cell_types::Union{Nothing, Vector{String}}
    domains::Union{Nothing, Vector{String}}
    coordinates::Vector{Tuple{Float64,Float64,Float64}}
    cell_distances::SparseMatrixCSC{Float64,Int64}
    max_diffusion::Float64 = 0.5
    density_damp::Float64 = 0.1
    diffusion_damp::Float64 = 0.3
    static::Bool = true
end
function Base.show(io::IO, ::MIME"text/plain",x::TissueState)
    println(io,"TissueState Struct:")
end

#-----------------------------------------------------------------------------#
# Temporal Struct - state progression info? keept track of part states?
# removing the saving part for now it is too heavy and not useful 
# during early proto
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct TemporalState
    total_steps::Union{Int64, Nothing} = nothing
    current_step::Int64 = 0
    sample_at::Union{Vector{Int64},Nothing} = nothing
end
function Base.show(io::IO, ::MIME"text/plain",x::TemporalState)
    println(io,"TemporalState Struct:")
end
#-----------------------------------------------------------------------------#
# Sample state - essentially collect all info to return a usable Sample
# also will contain initial condition information
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct SampleState
    n_cells::Int64 = 5000 
    n_genes::Int64 = 2000
    batch::Int64 = 1
    tissue::Union{TissueState, Nothing} = nothing
    cells::Union{Vector{CellState},Nothing} = nothing
    grn_set::Union{Dict, Nothing} = nothing
    gene_set::Union{GeneState, Nothing} = nothing
    temporal_state::Union{Nothing,TemporalState} = nothing
    biological_out::Vector{String} = ["rna"]
    out::Union{SparseMatrixCSC{Float64},Nothing} = nothing
end
function Base.show(io::IO, ::MIME"text/plain",x::SampleState)
    println(io,"SampleState Struct:")
end

#-----------------------------------------------------------------------------#
# Atlas state - creating collections of samples at the atlas level
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct AtlasState
    cells_per_sample::Tuple{Int64,Int64} = (1000, 10_000)
    genes_per_sample::Tuple{Int64,Int64} = (100, 10_000)
    n_samples::Int64 = 1
    n_batches::Int64 = 1
    domain_types::Union{Tuple{String, Int64, Int64},Vector{Tuple{String, Int64, Int64}}} = ("circle", 5,0)
    intial_cell_type::Int64 = 10
    grn_type::Union{Vector{String},Vector{Matrix}} = ["reprisilator","template_randomized"]
    biological_out::Vector{String} = Vector("rna",1)
    mosaic::Bool = false
    biological_noise::Int64 = 10
    batch_noise::Int64 = 10
    technical_noise::Int64 = 10
    samples:: Vector{SampleState}
    save_sample::Bool = true
end

function Base.show(io::IO, ::MIME"text/plain",x::AtlasState)
    println(io,"AtlasState Struct:")
end


#-----------------------------------------------------------------------------#
# Biological reference
# Will use this for now - not ideal 
# Mainly the labels - there could be a more compact way of parsing more
# labels and using them dynamically.
#-----------------------------------------------------------------------------#
Base.@kwdef struct BioRef
    coordinates::Tuple{Float64,Float64,Float64}
    type::String = "rna"
    biological_measurement::Union{Matrix{Float64},SparseMatrixCSC{Float64}}
    cell_labels::Union{Union{Vector{Int64}, Vector{String}}, Nothing}
    domain_labels::Union{Union{Vector{Int64}, Vector{String}}, Nothing}
    condition_labels::Union{Union{Vector{Int64}, Vector{String}}, Nothing}
end

function Base.show(io::IO, ::MIME"text/plain",x::BioRef)
    println(io,"BioRef Struct:")
end