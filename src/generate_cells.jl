

function add_shift_set(
    types::Vector{Int64},
    base_cell::BaseCells,
    shift_type::Symbol)::Dict
    shift_set = Dict{eltype(types), Tuple{NamedTuple, Vector{Int64}, Vector{Int64}}}()
    params = getfield(base_cell, shift_type)
    n_genes = base_cell.n_genes
    for t in types
        local_set = rand(1:n_genes, params.set_size)
        local_shift = add_shift(local_set, params.max_shift, n_genes)
        shift_set[t] = (params,local_set,local_shift)
    end
    return shift_set
end

# function add_chain_set(
#     types::Vector{Int64},
#     base_cell::BaseCells,
#     shift_type::Symbol)::Dict
#     shift_set = Dict{eltype(types), Dict{Int64,Tuple{NamedTuple, Vector, Vector}}}()
#     params = getfield(base_cell, shift_type)
#     n_genes = base_cell.n_genes
#     chain = 
#     for t in type

#     end
# end

function add_shift(
    gene_set::Vector{Int64},
    max_shift::Int64 = 500,
    n_genes::Int64 = 2000)::Vector{Int64}
    shifts = Vector{Int64}(undef,length(gene_set))
    for g in eachindex(gene_set)
        l_shift = rand(-max_shift:max_shift)
        if gene_set[g] + l_shift < 1
            l_shift = -gene_set[g] + 1
        end
        if gene_set[g] + l_shift > n_genes
            l_shift =  n_genes - gene_set[g]
        end
        shifts[g] = l_shift
    end
    return shifts
end

function add_cells(
    cell_params::BaseCells,
    tissue::Tissue;
    shift_type::Vector{String} = ["type_shifts","contact_shifts","wave_shifts"]
    )::BaseCells
    if shift_type isa String
        shift_type = [shift_types]
    end
    types = sort(unique(collect(values(tissue.cell_labels))))
    for shift in shift_type
        shift_set = add_shift_set(types, cell_params, Symbol(shift))
        setfield!(cell_params, Symbol(shift), shift_set)
    end
    cell_params.types = types
    cell_params.shift_type = shift_type
    return cell_params
end

#-----------------------------------------------------------------------------#
# New cells struct usage 
# Each cell is an individual agent
#-----------------------------------------------------------------------------#

function initialize_cells(
    tissue::TissueState,
    n_genes::Int64 = 2000,
    n_cells::Int64 = 5000)::Vector{CellState}
    #-------------------------------------------------------------------------#
    # pull and check
    #-------------------------------------------------------------------------#
    coordinates = tissue.coordinates
    if isnothing(tissue)
        throw(ArgumentError("Sample does not contain any tissue \n 
        Please run initialize_tissue first"))
    end
    cell_vec = Vector{CellState}(undef,n_cells)
    for c in 1:n_cells
        cell = initialize_cell(
            coordinates[c],
            n_genes)
        cell_vec[c] = cell
    end
    return cell_vec
end

function initialize_cell(
    coordinates::Tuple{Float64,Float64,Float64},
    n_genes::Int64 = 2000)::CellState
    rna = initialize_rank(n_genes)
    protein = initialize_rank(n_genes)
    cycle_position = rand(temporal_state) # 6 layers at the moment
    cell_info = CellInfo()
    #-------------------------------------------------------------------------#
    # build cell state struct
    #-------------------------------------------------------------------------#
    cell = CellState(
        cell_info = cell_info,
        cycle_position = cycle_position,
        ecosystem = nothing,
        coordinates = coordinates,
        chromatin_state = nothing,
        binding_state = nothing,
        rna_state = rna,
        protein_state = protein,
        metabolome_state = nothing,
        messaging_state = nothing)
    return cell
end

function add_cells(
    sample::SampleState,
    n_types::Int64 = 5,
    domain_bound::Bool = false)::SampleState
    if domain_bound
        
    else

    end
end


function add_cells_old(sample::SampleState)
    #-------------------------------------------------------------------------#
    # Initialize state vectors
    # Temporal state represent which biological level should we start the loop 
    # i.e. chromatin, tf_binding, etc. 
    # We will set a convention that state 1 is RNA
    #-------------------------------------------------------------------------#
    cell_info = [cell_type, domain]
    grn_local = Dict(k=> grn_set[k] for k in cell_info if haskey(grn_set, k))
    chromatin = initialize_state(grn_local,
        n_genes,
        :chromatin_remodelling,
        x -> x != 0)
    tf = initialize_state(grn_local,
        n_genes,
        :tf_binding,
        x -> x > 0)
    
    messaging = initialize_state(grn_local,
        n_genes,
        :messaging_output,
        x -> x != 0)
    metabolome = initialize_state(grn_local,
        n_genes,
        :metabolic_output,
        x -> x != 0)
end


using SparseArrays
function initialize_state(grn_set::Dict,
    n_genes::Int64,
    layer::Symbol,
    condition::Function)
    locs = Vector{Vector}(undef, length(grn_set))
    values = Vector{Vector}(undef, length(grn_set))
    count = 1
    for (_,grn) in grn_set
        l, v = grn_search(grn, layer, condition)
        locs[count] = l
        values[count] = v
        count += 1
    end
    concat_locs = vcat(locs)
    sparse_locs = unique(i -> concat_locs[i], eachindex(concat_locs))
    sparse_values = vcat(values)[sparse_locs]
    state = sparsevec(sparse_locs,sparse_values,n_genes)
    return state
end


function initialize_rank(n_genes)
    # for now we will keep this simple and not use sparse arrays
    return 1:n_genes
end


