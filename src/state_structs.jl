using SparseArrays
#-----------------------------------------------------------------------------#
# Gene struct - no need for mutability once set. 
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct GeneState
    n_genes::Int64 = 2000
    genes::Union{Vector{int64}, Vector{Str}}
    saturation_rank::Vector{Int64}
    leaky_rank::Vector{Int64}
    leak_rate::Vector{Int64} # if 0 then no change in baseline rank if more then it will go up till saturation of not repressed
    decay_rate::Vector{Float64}
    translation_efficiency::Vector{Float64}
    grn_state::Matrix{Int64}
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
    cellular_output::Vector{Int64}# Messages that will be diffused out
    metabolic_output::Vector{Int64}# Catch all for everything else
    chromatin_remodelling::Vector{Int64} # If they remodel chromatin - where do they do that?
    tf_binding::Vector{Int64}# If they bind somewhere where do they do that?
end

function Base.show(io::IO, ::MIME"text/plain",x::GRN)
    println(io,"GRN Struct")
end

#-----------------------------------------------------------------------------#
# Cell States struct
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct CellState
    temporal_state::Int64
    ecosystem::Matrix{Foat64}
    coordinates::Tuple{Float64,Float64,Float64}
    chromatin_state::Vector{Float64}
    binding_state::Vector{Float64}
    rna_state::Vector{Int64}
    protein_state::Vector{Int64}
    meta_state::Dict
    wrap_state::Dict
    trajectory::Vector{Str}
    grn_state::Dict
    cell_type::Union{Int32, Str}
    domain::Union{Int32, Str}
end

function Base.show(io::IO, ::MIME"text/plain",x::CellState)
    println(io,"CellState Struct:")
end

#-----------------------------------------------------------------------------#
# Tissue struct - general tissue information 
# Starting state I guesss since we want to add dynamic cells and mvt
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct TissueState
    n_cells::Int64 = 5000
    coordinates::Vector{Tuple{Float64,Float64,Float64}}
    cell_types::Vector{Int32}
    territories::Vector{Int32}
    cell_distance::SparseMatrixCSC{Float64}
    max_diffusion::Float64
    density_damp::Float64 = 0.1
    diffusion_damp::Float64 = 0.3
    static::Bool = true
end
function Base.show(io::IO, ::MIME"text/plain",x::TissueState)
    println(io,"TissueState Struct:")
end

#-----------------------------------------------------------------------------#
# Temporal Struct - state progression info? keept track of part states?
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct TemporalState
    total_steps::Int64
    current_step::Int64
    sample_at::Vector{Int64}
    past_state::Tuple{Vector{CellState},Vector{TissueState}}
    current_state::Tuple{Vector{CellState},Vector{TissueState}}
    save_states::Bool = false
end
function Base.show(io::IO, ::MIME"text/plain",x::TemporalState)
    println(io,"TemporalState Struct:")
end
#-----------------------------------------------------------------------------#
# Sample state - essentially collect all info to return a usable Sample
# also will contain initial condition information
#-----------------------------------------------------------------------------#
Base.@kwdef struct SampleState
    n_cells::Int64 = 5000 
    n_genes::Int64 = 2000
    batch::Int64
    tissue::TissueState
    cells::Vector{CellState}
    temporal_state::TemporalState
    spatial_linkage::Float64 = 0.3
    biological_out::Vector{Str} = Vector("rna",1)
    out::SparseMatrixCSC{Float64}
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
    grn_type::Union{Vector{str},Vector{Matrix}} = ["reprisilator","template_randomized"]
    biological_out::Vector{Str} = Vector("rna",1)
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
#-----------------------------------------------------------------------------#
Base.@kwdef struct BioRef
    coordinates::Tuple{Float64,Float64,Float64}
    type::Str = "rna"
    biological_measurement::Union{Matrix{Float64},SparseMatrixCSC{Float64}}
end

function Base.show(io::IO, ::MIME"text/plain",x::BioRef)
    println(io,"BioRef Struct:")
end