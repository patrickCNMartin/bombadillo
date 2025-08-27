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
        if any(x-> x > 1 || x < -1,sparse_values)
            print(sparse_values)
        end
        init_state = sparsevec(sparse_locs, sparse_values, n_genes)
    end
    return init_state
end


function cycle_chromatin!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)[1]
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
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)[1]
    chromatin_state = cell.chromatin_state
    tf_binding = cell.binding_state
    idx , val = findnz(tf_binding)
    for (i, v) in zip(idx, val)
        cs_prob = chromatin_state[i]
        tf_prob = abs(v) # we use abs since the signs only define if the gene is repressed or activated
        prob = (cs_prob + tf_prob) / 2  
        tf_binding[i] = sign(v) * prob
    end
    cell.binding_state = tf_binding
    return cell
end
function cycle_rna!(cell::CellState, gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)[1]
    rna = cell.rna_state
    saturation = gene_state.saturation_rank
    leak_rate = gene_state.leak_rate
    decay_rate = gene_state.decay_rate
    tf_binding = cell.binding_state
    n_genes = gene_state.n_genes
    idx, val = findnz(tf_binding)
    transcribe = [true, false]
    for (i, v) in zip(idx, val)
        # dirty trick if it comes at 0 just add one rank
        shift_max = abs(saturation[i] - rna[i]) + 1
        # I don't think random shift is the way to go here
        # the other issue is that this will make It
        # harder when adding decaying shifts
        # Will need to rework that 
        # TODO rework this while section 
        # TODO. add seperate functions for shifting
        # I think using (-) signs here makes more sense
        # It's more intuitive to think of decay as going down
        # even though it terms of rank index that means going up
        shift = -(sign(v) * round(logistic_sampling((0.0,Float64(shift_max))))) -
            leak_rate[i] -
            decay_rate[i]
        #shift =  -(sign(v) * shift_max) - leak_rate[i] - decay_rate[i]
        shift_prob = [abs(v), 1 - abs(v)] 
        shift_dist = Categorical(shift_prob)
        bound_site = transcribe[rand(shift_dist)]
        if bound_site
            rank = rna[i] - shift
            if 1 <= rank <= n_genes
                rank = maximum([saturation[i], rank])
            elseif rank < 1
                rank = 1
            else
                rank = n_genes
            end
            rna[i] = rank
        end
        
    end
    cell.rna_state = rna
    return cell
end
function cycle_protein!(cell::CellState,
    gene_state::GeneState)::CellState
    cell.cycle_position = circshift(temporal_state, (-1) * cell.cycle_position)[1]
    protein = cell.protein_state
    rna = cell.rna_state
    decay_rate = gene_state.decay_rate
    translation_efficiency = gene_state.translation_efficiency
    for p in eachindex(protein)
        rank = round(Int,rna[p] * (1 / translation_efficiency[p])) - decay_rate[p]
        protein[p] = rank
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
    # This will need to be update - using a function once we have 
    # more complex GRNs and GRN shifts
    info = [cell.cell_info.celltype,cell.cell_info.domain]
    grn_local = Dict(k => grn_set[k] for k in info if haskey(grn_set,k))
    saturation = gene_state.saturation_rank
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
    reg_index = abs.(regulation[:,2]) .∈ abs.(grn_layer)
    regulation = regulation[reg_index,:]

    for reg in eachrow(regulation)
        regulator = reg[2]
        target = reg[3]
        p_direction = sign(regulator)
        rank = protein[regulator]
        local_saturation = saturation[regulator]
        sat = (n_genes - rank) / (n_genes - local_saturation[p])
        if sat > 1.0 || sat < -1
                sat = 1.0
        elseif sat == 0.0
            # just add some noise here to avoid 0 sampling
            # Arbitratry but small chance of having spontanous activation
            sat = 0.05
        end
        prob = logistic_sampling((0.0, abs(sat)))
        cell_state[target] = p_direction * prob
    end
    setfield!(cell, cell_layer, cell_state)
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
    grn_set = sample.grn_set
    for cell in cells 
        state = cell.cycle_position
        cycle_cell[state](cell, gene_state)
        cycle_regulation!(cell,gene_state,grn_set)
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


function pull_state(cell::CellState, layer::Symbol, index::Int)
    return getfield(cell,layer)[index]
end