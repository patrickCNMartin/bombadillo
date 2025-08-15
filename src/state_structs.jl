#-----------------------------------------------------------------------------#
# Define const
#-----------------------------------------------------------------------------#
const temporal_state = [1,2,3,4,5,6]
#-----------------------------------------------------------------------------#
# depends
#-----------------------------------------------------------------------------#
using SparseArrays
using Parameters
using Printf
#-----------------------------------------------------------------------------#
# Gene struct - mutability for cases were the GRN types need to overwritw
# the gene params - eg high leak rate for repressiilator
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct GeneState
    n_genes::Int64 = 2000
    genes::Union{Vector{Int64}, Vector{String}}
    saturation_rank::Vector{Int64}
    leak_rate::Vector{Int64} # if 0 then no change in baseline rank if more then it will go up till saturation of not repressed
    decay_rate::Vector{Int64}
    translation_efficiency::Vector{Float64}
    regulator_strength::Vector{Float64}
end
function Base.show(io::IO, ::MIME"text/plain", x::GeneState)
    println(io, "GeneState Struct")
    println(io, "Number of Genes: ", x.n_genes)

    max_display = min(5, x.n_genes)  # Show at most 5 rows
    show_partial = x.n_genes > max_display

    headers = [
        "Gene",
        "Saturation Rank",
        "Leak Rate",
        "Decay Rate",
        "Translation Efficiency",
        "Regulator Strength"
    ]

    # Helper to format values based on type
    format_val(v) = v isa AbstractFloat ? @sprintf("%.3f", v) : string(v)

    # Extract rows with auto-formatting
    rows = [
        (
            format_val(x.genes[i]),
            format_val(x.saturation_rank[i]),
            format_val(x.leak_rate[i]),
            format_val(x.decay_rate[i]),
            format_val(x.translation_efficiency[i]),
            format_val(x.regulator_strength[i])
        )
        for i in 1:max_display
    ]

    # Compute column widths
    col_widths = [maximum(length, [h; getindex.(rows, j)]) for (j, h) in enumerate(headers)]

    # Print header
    for (j, h) in enumerate(headers)
        print(io, rpad(h, col_widths[j] + 2))
    end
    println(io)

    # Print rows
    for row in rows
        for (j, val) in enumerate(row)
            print(io, rpad(val, col_widths[j] + 2))
        end
        println(io)
    end

    # Indicate if truncated
    if show_partial
        println(io, "...")
    end
end

#-----------------------------------------------------------------------------#
# GRN struct - we use this to define GRN types - complex regulatory patterns
# E.g how do we encode a reprisilator? Will need to modify this
# No need for mutability once set.
#-----------------------------------------------------------------------------#
Base.@kwdef struct GRN
    regulatory_rel::Matrix{Int64}# columns as GRN "steps"
    regulators::Vector{Int64} # What genes are involved?
    strength_range::Tuple{Float64,Float64} = (0.5,1.0)
    regulator_strength::Vector{Float64} # How strongly are they "on" - master regulator always 1
    messaging_output::Union{Nothing,Vector{Int64}} = nothing# Messages that will be diffused out
    metabolic_output::Union{Nothing,Vector{Int64}} = nothing# Catch all for everything else
    chromatin_remodelling::Union{Nothing,Vector{Int64}} = nothing # If they remodel chromatin - where do they do that?
    tf_binding::Union{Nothing,Vector{Int64}} = nothing# If they bind somewhere where do they do that?
end

function Base.show(io::IO, ::MIME"text/plain", x::GRN)
    println(io, "GRN Struct")

    max_display = 5

    # Helper to format values
    format_val(v) = v isa AbstractFloat ? @sprintf("%.3f", v) :
                    v isa Nothing ? "—" : string(v)

    # Show regulatory_rel (matrix)
    n_rows, n_cols = size(x.regulatory_rel)
    println(io, "Regulatory Relations: $(n_rows)x$(n_cols) matrix")
    max_r = min(max_display, n_rows)
    max_c = min(max_display, n_cols)
    for i in 1:max_r
        row_str = join(format_val.(x.regulatory_rel[i, 1:max_c]), "  ")
        if n_cols > max_c
            row_str *= "  ..."
        end
        println(io, "  ", row_str)
    end
    if n_rows > max_r
        println(io, "  ...")
    end

    # Regulators
    println(io, "Regulators: ",
        x.regulators[1:min(max_display, length(x.regulators))],
        length(x.regulators) > max_display ? " ..." : "")

    # Strength range
    println(io, "Strength Range: ", "(", format_val(x.strength_range[1]),
        ", ", format_val(x.strength_range[2]), ")")

    # Regulator strength
    println(io, "Regulator Strength: ",
        [format_val(v) for v in x.regulator_strength[1:min(max_display, length(x.regulator_strength))]],
        length(x.regulator_strength) > max_display ? " ..." : "")

    # Optional vector fields
    for (label, vec) in [
        ("Messaging Output", x.messaging_output),
        ("Metabolic Output", x.metabolic_output),
        ("Chromatin Remodelling", x.chromatin_remodelling),
        ("TF Binding", x.tf_binding)
    ]
        if vec === nothing
            println(io, label, ": —")
        else
            println(io, label, ": ",
                vec[1:min(max_display, length(vec))],
                length(vec) > max_display ? " ..." : "")
        end
    end
end

#-----------------------------------------------------------------------------#
# CellInfo - meta data storage for fast lookup and replacement
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct CellInfo
    celltype::Union{String, Nothing} = nothing
    domain::Union{String, Nothing} = nothing
end

function Base.show(io::IO, ::MIME"text/plain", x::CellInfo)
    println(io, "CellInfo Struct")
    println(io, "Cell Type: ", x.celltype === nothing ? "—" : x.celltype)
    println(io, "Domain: ", x.domain === nothing ? "—" : x.domain)
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

function Base.show(io::IO, ::MIME"text/plain", x::CellState)
    println(io, "CellState Struct")

    # Cell Info
    if x.cell_info === nothing
        println(io, "Cell Info: —")
    else
        show(io, MIME"text/plain"(), x.cell_info)
    end

    # Cycle position
    println(io, "Cycle Position: ", x.cycle_position)

    # Ecosystem matrix
    if x.ecosystem === nothing
        println(io, "Ecosystem: —")
    else
        n_rows, n_cols = size(x.ecosystem)
        println(io, "Ecosystem: $(n_rows)x$(n_cols) matrix (showing up to 5x5)")
        max_r = min(5, n_rows)
        max_c = min(5, n_cols)
        for i in 1:max_r
            row_str = join((@sprintf("%.3f", x.ecosystem[i, j]) for j in 1:max_c), "  ")
            if n_cols > max_c
                row_str *= "  ..."
            end
            println(io, "  ", row_str)
        end
        if n_rows > max_r
            println(io, "  ...")
        end
    end

    # Coordinates
    println(io, "Coordinates: ", @sprintf("(%.3f, %.3f, %.3f)",
        x.coordinates[1], x.coordinates[2], x.coordinates[3]))

    # Helper for sparse/vector states
    function print_state(label, state)
        if state === nothing
            println(io, label, ": —")
        elseif state isa SparseVector
            println(io, label, ": SparseVector(len=", length(state),
                    ", nnz=", nnz(state), ")")
        elseif state isa Vector
            max_display = min(5, length(state))
            vals = state[1:max_display]
            println(io, label, ": ", vals, length(state) > max_display ? " ..." : "")
        else
            println(io, label, ": ", state)
        end
    end

    print_state("Chromatin State", x.chromatin_state)
    print_state("Binding State", x.binding_state)
    print_state("RNA State", x.rna_state)
    print_state("Protein State", x.protein_state)
    print_state("Metabolome State", x.metabolome_state)
    print_state("Messaging State", x.messaging_state)
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
function Base.show(io::IO, ::MIME"text/plain", x::TissueState)
    println(io, "TissueState Struct")

    # Cell types
    if x.cell_types === nothing
        println(io, "Cell Types: —")
    else
        max_display = min(5, length(x.cell_types))
        println(io, "Cell Types: ", x.cell_types[1:max_display],
            length(x.cell_types) > max_display ? " ..." : "")
    end

    # Domains
    if x.domains === nothing
        println(io, "Domains: —")
    else
        max_display = min(5, length(x.domains))
        println(io, "Domains: ", x.domains[1:max_display],
            length(x.domains) > max_display ? " ..." : "")
    end

    # Coordinates
    max_display = min(5, length(x.coordinates))
    coords_str = [@sprintf("(%.3f, %.3f, %.3f)", c[1], c[2], c[3]) for c in x.coordinates[1:max_display]]
    println(io, "Coordinates: ", coords_str,
        length(x.coordinates) > max_display ? " ..." : "")

    # Cell distances
    println(io, "Cell Distances: SparseMatrix($(size(x.cell_distances,1))x$(size(x.cell_distances,2))), nnz=", nnz(x.cell_distances))

    # Scalars
    println(io, "Max Diffusion: ", @sprintf("%.3f", x.max_diffusion))
    println(io, "Density Damp: ", @sprintf("%.3f", x.density_damp))
    println(io, "Diffusion Damp: ", @sprintf("%.3f", x.diffusion_damp))
    println(io, "Static: ", x.static)
end


#-----------------------------------------------------------------------------#
# Temporal Struct - state progression info? keept track of part states?
# removing the saving part for now it is too heavy and not useful 
# during early proto
#-----------------------------------------------------------------------------#
Base.@kwdef mutable struct TemporalState
    total_steps::Union{Int64, Nothing} = 0
    current_step::Int64 = 0
    sample_at::Union{Vector{Int64},Nothing} = nothing
end
function Base.show(io::IO, ::MIME"text/plain", x::TemporalState)
    println(io, "TemporalState Struct")

    println(io, "Total Steps: ", x.total_steps === nothing ? "—" : x.total_steps)
    println(io, "Current Step: ", x.current_step)

    if x.sample_at === nothing
        println(io, "Sample At: —")
    else
        max_display = min(5, length(x.sample_at))
        println(io, "Sample At: ", x.sample_at[1:max_display],
            length(x.sample_at) > max_display ? " ..." : "")
    end
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
    gene_state::Union{GeneState, Nothing} = nothing
    temporal_state::Union{Nothing,TemporalState} = nothing
    biological_out::Vector{String} = ["rna"]
    out::Union{SparseMatrixCSC{Float64},Nothing} = nothing
end
function Base.show(io::IO, ::MIME"text/plain", x::SampleState)
    println(io, "SampleState Struct")

    # Basic counts
    println(io, "Number of Cells: ", x.n_cells)
    println(io, "Number of Genes: ", x.n_genes)
    println(io, "Batch: ", x.batch)

    # Tissue
    if x.tissue === nothing
        println(io, "Tissue: —")
    else
        println(io, "Tissue: (see TissueState below)")
        show(io, MIME"text/plain"(), x.tissue)
    end

    # Cells
    if x.cells === nothing
        println(io, "Cells: —")
    else
        println(io, "Cells: Vector{CellState} (n=", length(x.cells), ")")
        if !isempty(x.cells)
            println(io, "  First cell preview:")
            show(io, MIME"text/plain"(), x.cells[1])
        end
    end

    # GRN set
    if x.grn_set === nothing
        println(io, "GRN Set: —")
    else
        println(io, "GRN Set: Dict with ", length(x.grn_set), " entries")
    end

    # Gene state
    if x.gene_state === nothing
        println(io, "Gene State: —")
    else
        println(io, "Gene State: (see GeneState below)")
        show(io, MIME"text/plain"(), x.gene_state)
    end

    # Temporal state
    if x.temporal_state === nothing
        println(io, "Temporal State: —")
    else
        println(io, "Temporal State: (see TemporalState below)")
        show(io, MIME"text/plain"(), x.temporal_state)
    end

    # Biological output
    max_display = min(5, length(x.biological_out))
    println(io, "Biological Output: ", x.biological_out[1:max_display],
        length(x.biological_out) > max_display ? " ..." : "")

    # Out matrix
    if x.out === nothing
        println(io, "Out Matrix: —")
    else
        println(io, "Out Matrix: SparseMatrix($(size(x.out,1))x$(size(x.out,2))), nnz=", nnz(x.out))
    end
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
function Base.show(io::IO, ::MIME"text/plain", x::AtlasState)
    println(io, "AtlasState Struct")

    # Basic ranges and counts
    println(io, "Cells per Sample: ", x.cells_per_sample)
    println(io, "Genes per Sample: ", x.genes_per_sample)
    println(io, "Number of Samples: ", x.n_samples)
    println(io, "Number of Batches: ", x.n_batches)

    # Domain types
    if x.domain_types isa Tuple
        println(io, "Domain Types: ", x.domain_types)
    else
        max_display = min(5, length(x.domain_types))
        println(io, "Domain Types: ", x.domain_types[1:max_display],
            length(x.domain_types) > max_display ? " ..." : "")
    end

    # Initial cell type
    println(io, "Initial Cell Type: ", x.intial_cell_type)

    # GRN type
    max_display = min(5, length(x.grn_type))
    println(io, "GRN Type: ", x.grn_type[1:max_display],
        length(x.grn_type) > max_display ? " ..." : "")

    # Biological output
    max_display = min(5, length(x.biological_out))
    println(io, "Biological Output: ", x.biological_out[1:max_display],
        length(x.biological_out) > max_display ? " ..." : "")

    # Booleans and noise
    println(io, "Mosaic: ", x.mosaic)
    println(io, "Biological Noise: ", x.biological_noise)
    println(io, "Batch Noise: ", x.batch_noise)
    println(io, "Technical Noise: ", x.technical_noise)

    # Samples
    println(io, "Samples: Vector{SampleState} (n=", length(x.samples), ")")
    if !isempty(x.samples)
        println(io, "  First sample preview:")
        show(io, MIME"text/plain"(), x.samples[1])
    end

    # Save sample flag
    println(io, "Save Sample: ", x.save_sample)
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

function Base.show(io::IO, ::MIME"text/plain", x::BioRef)
    println(io, "BioRef Struct")

    # Coordinates
    println(io, "Coordinates: ", @sprintf("(%.3f, %.3f, %.3f)", x.coordinates...))

    # Type
    println(io, "Type: ", x.type)

    # Biological measurement matrix
    n_rows, n_cols = size(x.biological_measurement)
    matrix_type = x.biological_measurement isa SparseMatrixCSC ? "SparseMatrix" : "Matrix"
    println(io, "Biological Measurement: ", matrix_type, " ($(n_rows)x$(n_cols)), nnz=", nnz(x.biological_measurement))

    # Helper for labels
    function print_labels(label, vec)
        if vec === nothing
            println(io, label, ": —")
        else
            max_display = min(5, length(vec))
            println(io, label, ": ", vec[1:max_display],
                length(vec) > max_display ? " ..." : "")
        end
    end

    print_labels("Cell Labels", x.cell_labels)
    print_labels("Domain Labels", x.domain_labels)
    print_labels("Condition Labels", x.condition_labels)
end
