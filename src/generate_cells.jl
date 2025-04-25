
function add_gene_set(
    types::AbstractVector;
    n_genes::Int64 = 2000,
    set_size::Int64 = 10)::Dict
    type_sets = Dict(k => Vector{Int}() for k in types)
    for t in types
        type_sets[t] = rand(1:n_genes, set_size)
    end
    return type_sets
end

function add_shift(
    gene_set::Dict;
    max_shift::Int64 = 500,
    n_genes::Int64 = 2000)::Dict
    shifts = Dict{Any, Vector{Int}}()
    for g in keys(gene_set)
        l_genes = collect(values(gene_set[g]))
        l_shifts = Vector{Int}(undef, length(l_genes))

        for v in eachindex(l_genes)
            l_shift = rand(-max_shift:max_shift)
            if l_genes[v] + l_shift < 1
                l_shift = -l_genes[v] + 1
            end
            if l_genes[v] + l_shift > n_genes
                l_shift =  n_genes - l_genes[v]
            end
            
            l_shifts[v] = l_shift
        end
        shifts[g] = l_shifts
    end

    return shifts
end

function add_cells(
    tissue::Tissue;
    n_genes::Int64 = 2000,
    set_size::Int64 = 10,
    max_shift::Int64 = 500,
    max_diffusion::Int64 = 5,
    dampening::Float64 = 0.1)::Tissue
    types = sort(unique(collect(values(tissue.cell_labels))))
    type_sets = add_gene_set(types, n_genes = n_genes, set_size = set_size)
    type_shifts = add_shift(type_sets, max_shift = max_shift, n_genes = n_genes)
    contact_sets = add_gene_set(types, n_genes = n_genes, set_size = set_size)
    contact_shifts = add_shift(contact_sets, max_shift = max_shift, n_genes = n_genes)
    wave_sets = add_gene_set(types, n_genes = n_genes, set_size = set_size)
    wave_shifts = add_shift(wave_sets, max_shift = max_shift, n_genes = n_genes)
    base_cell = BaseCellType(
        types = types,
        type_sets = type_sets,
        type_shifts = type_shifts,
        contact_sets = contact_sets,
        contact_shifts = contact_shifts,
        wave_sets = wave_sets,
        wave_shifts = wave_shifts,
        wave_diffusion = max_diffusion,
        wave_damp = dampening)
    tissue.base_cell_types = base_cell
    return tissue
end

