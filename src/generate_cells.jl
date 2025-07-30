#-----------------------------------------------------------------------------#
# Setting a constant for temporal state 
# we might reuse this at different places and it would be cleaner to change
# it here only.
#-----------------------------------------------------------------------------#
const temporal_state = [1,2,3,4,5,6]


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

function initialize_cell(
    cell_type::Union{Int64, String},
    domain::Union{Int64, String},
    grn_set::Dict,
    n_genes::Int64 = 2000,
    coordinate_range::Tuple{Float64,Float64} = (0.0,1.0))::CellState
    #-------------------------------------------------------------------------#
    # Initialize coordinates
    #-------------------------------------------------------------------------#
    x = generate_coordinates(coordinate_range)
    y = generate_coordinates(coordinate_range)
    z = generate_coordinates(coordinate_range)
    coordinates = tuple.(x,y,z)
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
    rna = initialize_rank(n_genes)
    protein = initialize_rank(n_genes)
    messaging = initialize_state(grn_local,
        n_genes,
        :messaging_output,
        x -> x != 0)
    metabolome = initialize_state(grn_local,
        n_genes,
        :metabolic_output,
        x -> x != 0)
    cycle_position = rand(temporal_state) # 6 layers at the moment 
    #-------------------------------------------------------------------------#
    # build cell state struct
    #-------------------------------------------------------------------------#
    cell = CellState(
        cell_type = cell_type,
        domain = domain,
        cycle_position = cycle_position,
        grn_set = grn_local,
        ecosystem = nothing,
        coordinates = coordinates,
        chromatin_state = chromatin,
        binding_state = tf,
        rna_state = rna,
        protein_state = protein,
        metabolome_state = metabolome,
        messaging_state = messaging)
    return cell
end


function generate_coordinates(coordinate_range::Tuple{Float64,Float64} = (0.0,1.0))
    return rand(Uniform(coordinate_range[1],coordinate_range[2]))
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


