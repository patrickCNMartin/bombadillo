using Statistics
function repressilator(
    regulator_strength::Vector{Float64},
    remodeler_strength::Vector{Float64},
    gene_state::GeneState;
    n_regulators::Int64 = 3,
    percentile::Float64 = 0.1)::Tuple{GRN,GeneState}
    #-------------------------------------------------------------------------#
    # Select some regulators 
    # One thing that is required is that we cannot randomly select 
    # Only one has to be on. The other turns off. The bahaviors needs to be 
    # forced. Always start with a single repressor on
    # Will need to throw in some check here later
    # Will do when we create the GRN creator tool kit
    #-------------------------------------------------------------------------#
    regulators = check_selected_regulators(regulator_strength)
    remodelers = check_selected_regulators(remodeler_strength)
    regulators = intersect(regulators,remodelers)
    expressed = quantile(regulators, percentile)
    expressed_regs = rand(regulators[regulators .< expressed], 1)
    repressed = quantile(regulators, 1 - percentile)
    repressed_regs = rand(regulators[regulators .> repressed], n_regulators - 1)
    # local_reg = regulators[(regulators .> n_regulators) .& (regulators .< gene_state.n_genes)]
    # loc = rand(local_reg)
    # all_regulators = local_reg[loc:(loc + (n_regulators - 1))]
    #-------------------------------------------------------------------------#
    # Function tp shift indices to make a circular repressilator of 
    # arbitrary size - sign just add a negavtive sign for "repression"
    # Returns Matrix{Int64}
    # Multiply by -1 to define that there is a repressive relationship
    #-------------------------------------------------------------------------#
    all_regulators = vcat(expressed_regs, repressed_regs)
    reg_rel = cyclic_permuations(all_regulators)
    reg_rel = reg_rel .* (-1)
    #-------------------------------------------------------------------------#
    # Define strength of regulators
    # Repressed ones are low
    # With negative strenghts to indicate repression
    # this could be redundant.
    #-------------------------------------------------------------------------#
    reg_strengths = (-1) .* vcat(1, rand(Uniform(0.0,percentile),n_regulators - 1))
    #reg_strengths = repeat([-1.0],n_regulators)
    remod_strength = repeat([1.0],n_regulators)
    #-------------------------------------------------------------------------#
    # Update gene state 
    #-------------------------------------------------------------------------#
    gene_state.regulator_strength[all_regulators] .= reverse(reg_strengths)
    gene_state.remodeler_strength[all_regulators] .= remod_strength
    #-------------------------------------------------------------------------#
    # For the repressillator there is no need to have any other output
    # Except for TF binding since they bind and inhibit expression.
    # Build final grn Struct
    # Currently unmatable struct since once the GRN is set the rest happens
    # elsewhere. This just serves a behavior template
    #-------------------------------------------------------------------------#
    grn = GRN(regulatory_rel = reg_rel,
        regulators = all_regulators,
        strength_range = (0.99,1.0),
        regulator_strength = reverse(reg_strengths),
        remodeler_strength = remod_strength,
        messaging_output = nothing,
        metabolic_output = nothing,
        chromatin_remodelling = all_regulators,
        tf_binding = (-1) * all_regulators)
    return grn, gene_state
end

function check_selected_regulators(
    regulator_strength::Vector{Float64})::Vector{Int64}
    return findall(x -> x == 0, regulator_strength)
end

## Not a great function either too complicated or too simple
function compute_grn_overlaps(
    regulator_strength::Vector{Float64};
    overlap_range::Tuple{Float64,Float64} = (0.0,0.0),
    strength_range::Tuple{Float64,Float64} = (0.5,1.0),
    g::Int64 = 3)::Tuple{Vector{Int64},Vector{Float64}}
    #-------------------------------------------------------------------------#
    # If we don't want overlaps we find that have strength of 0 
    # We assume that 0 is non activate genes so free to use
    #-------------------------------------------------------------------------#
    if overlap_range[1] == 0.0 && overlap_range[2] == 0.0
        locs = findall(regulator_strength .== 0.0)
        regulators = sample(locs, g, replace = false)
        strenghts = rand(Uniform(strength_range[1],strength_range[2]),g)
        regulator_strength[regulators] .= strenghts
    #-------------------------------------------------------------------------#
    # If we want overlaps we can define a overlap range 
    # From maybe you have an overlap to yes you definitely do have one 
    # The aim here is to hav GRN interactions
    # GRN are small entities that can have run-away and emmergent properties
    #-------------------------------------------------------------------------#
    else
        over_g = rand(Int(g * overlap_range[1]):Int(g * overlap_range[2]))
        new_g = g - over_g
        over_locs = findall(regulator_strength .> 0)
        new_locs = findall(regulator_strength .== 0)
        regulators = [sample(over_locs,over_g);sample(new_locs,new_g)]
        strenghts = rand(Uniform(strength_range[1],strength_range[2]), over_g + new_g)
        regulator_strength[regulators] .= strenghts
    end
    return regulators, strenghts
end


# function grn_search(grn::GRN,
#     field_name::Symbol,
#     condition::Function)::Tuple
#     field = getfield(grn, field_name)
#     locs = findall(condition, field)
#     values = field[locs]
#     return (locs, values)
# end

#-----------------------------------------------------------------------------#
# add grns for simple type such as cells or domain 
# we will add messaging_output seperately i guess
#-----------------------------------------------------------------------------#
function add_grns(sample::SampleState,
    use::Union{Symbol, Vector{Symbol}} = :cell_types,
    overwrite::Bool = false)::SampleState
    #-------------------------------------------------------------------------#
    # This will need to be refactored
    # It's not very elegant
    #-------------------------------------------------------------------------#
    if isa(use, Symbol)
        grns = unique(check_field_value(sample.tissue, use))
    else
        grns = unique([check_field_value(sample.tissue,i) for i in use])
    end
    #-------------------------------------------------------------------------#
    # Check if GRns are already present otherwise add to existing 
    #-------------------------------------------------------------------------#
    if isnothing(sample.grn_set)
        grn_set = Dict{String,GRN}()
    else
        grn_set = sample.grn_set
    end

    #-------------------------------------------------------------------------#
    # check if we overwrite or not - this is useful if we want to add
    # more than one grn per category
    #-------------------------------------------------------------------------#
    if !overwrite
        grns = make_unique_keys(grn_set, grns)
    end
    #-------------------------------------------------------------------------#
    # Finally we pull regulator strength - intilized at 0
    # The idea is that during initialization we just use the ones with 0s
    # we can add overlaps and the strength will determine how strongly
    # the regulator will resist change
    #-------------------------------------------------------------------------#
    gene_state = sample.gene_state
    regulator_strength = sample.gene_state.regulator_strength
    remodeler_strength = sample.gene_state.remodeler_strength
    for g in eachindex(grns)
        grn,gene_state = repressilator(
            regulator_strength,
            remodeler_strength,
            gene_state)
        grn_set[grns[g]] = grn
    end
    sample.grn_set = grn_set
    sample.gene_state = gene_state
    
    return sample
end


function add_grns!(sample::SampleState,
    use::Union{Symbol, Vector{Symbol}} = :cell_types,
    overwrite::Bool = false)::SampleState
    #-------------------------------------------------------------------------#
    # This will need to be refactored
    # It's not very elegant
    #-------------------------------------------------------------------------#
    if isa(use, Symbol)
        grns = unique(check_field_value(sample.tissue, use))
    else
        grns = unique([check_field_value(sample.tissue,i) for i in use])
    end
    #-------------------------------------------------------------------------#
    # Check if GRns are already present otherwise add to existing 
    #-------------------------------------------------------------------------#
    if isnothing(sample.grn_set)
        grn_set = Dict{String,GRN}()
    else
        grn_set = sample.grn_set
    end

    #-------------------------------------------------------------------------#
    # check if we overwrite or not - this is useful if we want to add
    # more than one grn per category
    #-------------------------------------------------------------------------#
    if !overwrite
        grns = make_unique_keys(grn_set, grns)
    end
    #-------------------------------------------------------------------------#
    # Finally we pull regulator strength - intilized at 0
    # The idea is that during initialization we just use the ones with 0s
    # we can add overlaps and the strength will determine how strongly
    # the regulator will resist change
    #-------------------------------------------------------------------------#
    gene_state = sample.gene_state
    regulator_strength = sample.gene_state.regulator_strength
    remodeler_strength = sample.gene_state.remodeler_strength
    for g in eachindex(grns)
        grn,gene_state = repressilator(
            regulator_strength,
            remodeler_strength,
            gene_state)
        grn_set[grns[g]] = grn
    end
    sample.grn_set = grn_set
    sample.gene_state = gene_state
    
    return sample
end
