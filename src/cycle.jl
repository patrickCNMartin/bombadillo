function let_live(
    sample::SampleState,
    cycles::Int = 100,
    sample_genes::Union{Nothing, String, Vector{String}} = nothing;
    pull_layer::Symbol = :rna_state
)
    
    if isnothing(sample.temporal_state)
        cycle_state = TemporalState(total_steps = cycles)
        initialize_cycling!(sample)
        sample.temporal_state = cycle_state
    end
    
    
    gene_cycle = initialize_gene_pull(sample_genes, cycles, sample.n_cells)
    
    
    for cycle in 1:cycles
        cycle_sample!(sample, cycle)
        println(string("Cycle:",cycle))
        if gene_cycle !== nothing
            gene_pull!(gene_cycle, sample_genes, sample, cycle, pull_layer)
        end
    end
    return isnothing(gene_cycle) ? sample : (sample, gene_cycle)
end


using DataFrames
function initialize_gene_pull(sample_genes, cycles::Int, n_cells::Int)::Union{Vector{DataFrame},Nothing}
    if isnothing(sample_genes)
        return nothing

    elseif isa(sample_genes, String)
        dfs = [DataFrame([Symbol("Cycle$(c)") => Vector{Float64}(undef, n_cells) 
                          for c in 1:cycles]...)]
        
    elseif isa(sample_genes, Vector{String})
        dfs = Vector{DataFrame}(undef, length(sample_genes))
        for (i, g) in enumerate(sample_genes)
            dfs[i] = DataFrame([Symbol("Cycle$(c)") => Vector{Float64}(undef, n_cells) 
                                for c in 1:cycles]...)
        end

    else
        throw(TypeError("Unknown type for calling genes - use String or Vector{String}"))
    end

    return dfs
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
        
        if field == :chromatin_remodelling
            str = grn.remodeler_strength
            reg = grn.regulators
        else
            str = grn.regulator_strength
            reg = grn.regulators
        end
        
        effect_strength = grn.strength_range
        layer = getfield(grn, field)
        # Don't initalize what is not present
        # Not all GRNs will have all fields present
        if isnothing(layer)
            continue
        end
        for pr in eachindex(layer)
            loc_key = abs(layer[pr])
            # Not sure if overwiting with stronger one is the way to go
            # maybe be mean strenght? 
            if haskey(state,loc_key) && loc_key ∉ reg
                st =  rand(Uniform(effect_strength[1],effect_strength[2]))
                if abs(state[loc_key]) < abs(st)
                    state[loc_key] =  sign(layer[pr]) * abs(st)
                end
            elseif haskey(state, loc_key) && loc_key ∈ reg
                st = only(str[reg .== loc_key])
                if abs(state[loc_key]) < abs(st)
                    state[loc_key] = sign(layer[pr]) * abs(st)
                end
            elseif !haskey(state,loc_key) && loc_key ∈ reg
                st = only(str[reg .== loc_key])
                state[loc_key] = st
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
        if any(x-> x > 1 || x < -1,sparse_values)
            print(sparse_values)
        end
        init_state = sparsevec(sparse_locs, sparse_values, n_genes)
    end
    return init_state
end


function cycle_chromatin_binding!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)[1]
    chromatin_state = cell.chromatin_state
    tf_binding = cell.binding_state
    idx, val = findnz(chromatin_state)
    for (i, v) in zip(idx,val)
        binding_prob = mean([abs(v),abs(tf_binding[i])]) * sign(tf_binding[i]) 
        tf_binding[i] = binding_prob
    end
    cell.binding_state = tf_binding
    return cell
end



using Random, Distributions
function cycle_rna!(cell::CellState, gene_state::GeneState)::CellState
    # advance cell cycle position
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)[1]

    # pull data
    expr = cell.rna_state                       
    leak_rate = gene_state.leak_rate            
    decay_rate = gene_state.decay_rate   
    tf_binding = cell.binding_state 
    n_genes = gene_state.n_genes          
    saturation = gene_state.saturation  

    # only consider nonzero TF bindings
    idx, val = findnz(tf_binding)
    transcribe_options = [true, false]

    @inbounds for (i, v) in zip(idx, val)
        # probabilistic transcription based on TF strength
        shift_prob = [abs(v), 1 - abs(v)]
        shift_dist = Categorical(shift_prob)
        bound_site = transcribe_options[rand(shift_dist)]
        local_decay = logistic_sampling((1.0, Float64(n_genes - saturation[i])),5.0) * decay_rate[i]
        local_decay = round(Int64, local_decay) + 1
        local_leak = logistic_sampling((1.0,Float64(n_genes - saturation[i])), 10.0) * leak_rate[i]
        local_leak = round(Int64, local_leak) + 1
        if bound_site
            v_loc = logistic_sampling((1.0, Float64(n_genes - saturation[i])),10.0)
            v_loc = round(Int64,  v_loc)
            delta = -local_leak + local_decay - v_loc 
        else
            delta = -local_leak + local_decay 
        end
        expr[i] += delta

        if saturation !== nothing
            expr[i] = min(expr[i], saturation[i])
        end
        expr[i] = max(expr[i], 0)
    end

    cell.rna_state = expr
    return cell
end


function cycle_protein!(cell::CellState,
    gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)[1]
    protein = cell.protein_state
    rna = cell.rna_state
    decay_rate = gene_state.decay_rate
    translation_efficiency = gene_state.translation_efficiency
    n_genes = gene_state.n_genes
    for p in eachindex(protein)
        prot_shift = n_genes - rna[p]
        decay_shift = round(Int64,prot_shift * decay_rate[p])
        translation_shift = round(Int64,prot_shift * translation_efficiency[p])
        protein_level = rna[p] + decay_shift + translation_shift
        if protein_level > n_genes
            protein_level = n_genes
        end
        protein[p] = protein_level
    end
    cell.protein_state = protein
    return cell
end



function cycle_messaging!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)[1]
    return cell
end

function cycle_metabolome!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)[1]
    return cell
end

function cycle_regulation!(cell::CellState,
    gene_state::GeneState,
    grn_set::Dict{String,GRN})::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)[1]
    # This will need to be update - using a function once we have 
    # more complex GRNs and GRN shifts
    info = [cell.cell_info.celltype,cell.cell_info.domain]
    grn_local = Dict(k => grn_set[k] for k in info if haskey(grn_set,k))
    saturation = gene_state.saturation
    n_genes = gene_state.n_genes
    protein = cell.protein_state
    for (_,g) in grn_local
        # assuming that there will only be TF binding and chromatin State
        if !isnothing(g.tf_binding)
            update_regulation!(cell,
                g,
                protein,
                saturation,
                n_genes,
                grn_layer = :tf_binding,
                cell_layer = :binding_state)
        end
        if !isnothing(g.chromatin_remodelling)
            update_regulation!(cell,
                g,
                protein,
                saturation,
                n_genes,
                grn_layer = :chromatin_remodelling,
                cell_layer = :chromatin_state)
        end
    end
    return cell
end

function update_regulation!(cell::CellState,
    grn::GRN,
    protein::Vector{Int64},
    saturation::Vector{Int64},
    n_genes::Int64;
    grn_layer::Symbol = :tf_binding,
    cell_layer::Symbol = :binding_state)::CellState
    grn_layer = getfield(grn, grn_layer)
    cell_state = getfield(cell, cell_layer)
    # pull GRN rels and subset 
    regulation = grn.regulatory_rel
    reg_index = indexin(abs.(regulation[:,2]),abs.(grn_layer))
    regulation = regulation[reg_index,:]
    
    for reg in eachrow(regulation)
        regulator = abs(reg[2])
        target = abs(reg[3])
        p_direction = sign(reg[2])
        rank = protein[regulator]
        local_saturation = saturation[regulator]
        max_saturation_delta = n_genes - local_saturation
        max_rank_delta = (n_genes - rank) + 1
        sat =  max_rank_delta / max_saturation_delta
        #println(rank - local_saturation)
        prob = p_direction * hill_signed(sat)
        #new_prob = clamp(cell_state[target] * prob, 0.001,0.999)
        new_prob = clamp(cell_state[target] + p_direction * prob, 0.1,0.999)
        cell_state[target] = new_prob  
    end
    setfield!(cell, cell_layer, cell_state)
    return cell
end


const cycle_cell = Dict{Int, Function}(
    1 => cycle_chromatin_binding!,
    2 => cycle_rna!,
    3 => cycle_protein!,
    4 => cycle_messaging!,
    5 => cycle_metabolome!)

function cycle_sample!(sample::SampleState, cycle::Int)::SampleState
    cells = sample.cells
    gene_state = sample.gene_state
    grn_set = sample.grn_set
    for cell in cells 
        state = cell.cycle_position
        if state != 6
            cycle_cell[state](cell, gene_state)
        else
            cycle_regulation!(cell,gene_state,grn_set)
        end 
    end
    sample.temporal_state.total_steps = cycle
    sample.temporal_state.current_step = cycle
    return sample
end



function gene_pull!(
    gene_cycle::Vector{DataFrame},
    genes::Union{String,Vector{String}},
    sample::SampleState,
    cycle::Int,
    pull_layer::Symbol)::Vector{DataFrame}

    gene_index = findall(g -> g ∈ genes, sample.gene_state.genes)
    cells = sample.cells

    for (i, gidx) in enumerate(gene_index)
        
        state = similar(cells, Float64)
        @inbounds @simd for j in eachindex(cells)
            state[j] = pull_state(cells[j], pull_layer, gidx)
        end
        @views gene_cycle[i][!, cycle] .= state
    end

    return gene_cycle
end


function pull_state(
    cell::CellState,
    layer::Symbol,
    index::Int)
    return getfield(cell,layer)[index]
end


# Moved from visualization.jl (no plotting dependency)
function delta_shift(pull::DataFrame)::Matrix
    n_cells = size(pull, 1)
    num_time = size(pull, 2)
    delta_mat = zeros(Float64, n_cells, num_time)
    for t in 2:num_time
        @views delta_mat[:, t] .= pull[:, t] .- pull[:, t-1]
    end
    return delta_mat
end
