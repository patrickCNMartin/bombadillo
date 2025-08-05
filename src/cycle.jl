function cycle_sample(
    sample::SampleState,
    cycles::Int64 = 100)::SampleState
    if isnothing(sample.temporal_state)
        temporal_state = TemporalState(
            total_steps = cycles)
        initialize_cycling!(sample)
    end
    for cycle in 1:cycles
        
    end

end


function initialize_cycling!(sample::SampleState)::SampleState
    cells = check_field_value(sample, :cells)
    cell_info = check_field_value.(cells, :cell_info)
    grn_set = check_field_value(sample,:grn_set)
    n_genes = sample.n_genes
    for cell in eachindex(cells)
        info = [cell_info[cell].celltype,cell_info[cell].domain]
        grn_local = Dict(k => grn_set[k] for k in info if haskey(grn_set,k))
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
        cells[cell].chromatin_state = chromatin
        cells[cell].binding_state = tf
        cells[cell].messaging_state = messaging
        cells[cell].metabolome_state = metabolome
    end
    sample.cells = cells
    return sample
end



using SparseArrays
function initialize_state(grn_set::Dict,
    n_genes::Int64,
    field::Symbol)
    state = Dict{Int64, Int64}()
    for (_, grn) in grn_set
        reg = grn.regulators
        str = grn.regulator_strength
        layer = getfield(grn, field)
        if isnothing(layer)
            contiue
        end
        
        for (r, s, l) in zip(str,layer)
            if haskey(state, reg)
                
            end
            index_to_value[idx] = val
        end
    end
    sparse_locs = collect(keys(index_to_value))
    sparse_values = collect(values(index_to_value))
    return sparsevec(sparse_locs, sparse_values, n_genes)::SparseVector{Int64, Int64}
end
