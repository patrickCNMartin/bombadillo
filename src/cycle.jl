function let_live(
    sample::SampleState,
    cycles::Int64 = 100)::SampleState
    if isnothing(sample.temporal_state)
        cycle_state = TemporalState(
            total_steps = cycles)
        initialize_cycling!(sample)
        sample.temporal_state = cycle_state
    end
    for cycle in 1:cycles
        cycle_sample!(sample, cycle)
    end
    return sample
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
            :chromatin_remodelling)
        tf = initialize_state(grn_local,
            n_genes,
            :tf_binding)
        messaging = initialize_state(grn_local,
            n_genes,
            :messaging_output)
        metabolome = initialize_state(grn_local,
            n_genes,
            :metabolic_output)
        cells[cell].chromatin_state = chromatin
        cells[cell].binding_state = tf
        cells[cell].messaging_state = messaging
        cells[cell].metabolome_state = metabolome
    end
    sample.cells = cells
    return sample
end



using SparseArrays
using Distributions
function initialize_state(grn_set::Dict{String, Bombadillo.GRN},
    n_genes::Int64,
    field::Symbol)::SparseVector{Float64, Int64}
    state = Dict{Int64, Float64}()
    for (_, grn) in grn_set
        reg = grn.regulators
        str = grn.regulator_strength
        effect_strength = grn.strength_range
        layer = getfield(grn, field)
        # Don't initalize what is not present
        # Not all GRNs will have all fields present
        if isnothing(layer)
            continue
        end
        for pr in eachindex(layer)
            loc_key = abs(layer[pr])
            if haskey(state,loc_key) && loc_key ∉ reg
                st =  rand(Uniform(effect_strength[1],effect_strength[2]))
                if state[loc_key] < st
                    state[loc_key] =  sign(layer[pr]) * st
                end
            elseif haskey(state, loc_key) && loc_key ∈ reg
                st = str[reg == loc_key]
                if state[loc_key] < st
                    state[loc_key] = sign(layer[pr]) * st
                end
            else
                state[loc_key] = sign(layer[pr]) * rand(Uniform(effect_strength[1],effect_strength[2]))
            end
        end
    end
    if isempty(state)
        init_state = spzeros(Float64, n_genes)
    else
        sparse_locs = collect(keys(state))
        sparse_values = collect(values(state))
        init_state = sparsevec(sparse_locs, sparse_values, n_genes)
    end
    return init_state
end


function cycle_chromatin!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)
    # what is the probablity of getting a state change
    chromatin_state = cell.chromatin_state
    idx, val = findnz(chromatin_state)
    flip = [1, -1]
    for (i, v) in zip(idx,val)
        prob = [abs(v), 1 - abs(v)]
        dist = Categorical(prob)
        chromatin_state[i] = flip[rand(dist)] * v
    end
    cell.chromatin_state = chromatin_state
    return cell
end

function cycle_binding!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)
    chromatin_state = cell.chromatin_state
    tf_binding = cell.binding_state
    idx , val = findnz(tf_binding)
    for (i, v) in zip(idx, val)
        cs_prob = chromatin_state[i]
        tf_prob = v
        prob = (cs_prob + tf_prob) / 2
        tf_binding[i] = sign(v) * prob
    end
    cell.binding_state = tf_binding
    return cell
end
function cycle_rna!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)
    rna = cell.rna_state
    saturation = gene_state.saturation_rank
    leak_rate = gene_state.leak_rate
    decay_rate = gene_state.decay_rate
    tf_binding = cell.binding_state
    idx, val = findnz(tf_binding)
    for (i, v) in zip(idx, val)
        # dirty trick if it comes at 0 just add one rank
        shift_max = abs(saturation[i] - rna[i]) + 1
        shift = (sign(v) * rand(0:shift_max)) +
            leak_rate[i] +
            decay_rate[i]
        rank = rna[i] + shift
        rna[i] = minimum([saturation[i], rank])
    end
    cell.rna_state = rna
    return cell
end
function cycle_protein!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)
    protein = cell.protein_state
    rna = cell.rna_state
    decay_rate = gene_state.decay_rate
    translation_efficiency = gene_state.translation_efficiency
    for p in eachindex(protein)
        rank = round(Int,rna[p] * translation_efficiency[p]) + decay_rate[p]
        protein[p] = rank
    end
    cell.protein_state = protein
    return cell
end
function cycle_messaging!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)
    return cell
end

function cycle_metabolome!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)
    return cell
end


const cycle_cell = Dict{Int, Function}(
    1 => cycle_chromatin!,
    2 => cycle_binding!,
    3 => cycle_rna!,
    4 => cycle_protein!,
    5 => cycle_messaging!,
    6 => cycle_metabolome!)

function cycle_sample!(sample::SampleState, cycle::Int)::SampleState
    cells = sample.cells
    gene_state = sample.gene_state
    for cell in cells 
        state = cell.cycle_position
        cycle_cell[state](cell, gene_state)
    end
    sample.temporal_state.total_steps = cycle
    sample.temporal_state.current_step = cycle
    return sample
end

