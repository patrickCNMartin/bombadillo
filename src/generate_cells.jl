
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
    cell_type::Union{Int32, String},
    domain::Union{Int32, String},
    grn_set::Vector{GRN},
    gene_state::GeneState,
    coordinate_range::Vector{Float64} = [0.0,1.0])::CellState
    #-------------------------------------------------------------------------#
    # First we generate coordinates
    #-------------------------------------------------------------------------#
    x = generate_coordinates(coordinate_range)
    y = generate_coordinates(coordinate_range)
    z = generate_coordinates(coordinate_range)
end


function generate_coordinates(coordinate_range::Vector{Float64} = [0.0,1.0])
    return rand(Uniform(coordinate_range[1],coordinate_range[2]))
end