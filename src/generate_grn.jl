
function repressilator(
    regulator_strength::Vector{Float64},
    gene_state::GeneState;
    n_regulators::Int64 = 3)::Tuple{GRN,GeneState}
    #-------------------------------------------------------------------------#
    # This function selects genes to use in the repressilator
    # We can define overlaps with other genes 
    # And how strong they are "on"
    # Return Vector{Int64} and Vector{Float64}
    # Hardcoded values for now just to test it
    # I think I would need a GRN designed model
    #-------------------------------------------------------------------------#
    regulators, strengths = compute_grn_overlaps(
        regulator_strength,
        overlap_range = (0.0,0.0),
        strength_range = (0.99,1.0),
        g = n_regulators)
    gene_state.leak_rate[regulators] .= 200
    #-------------------------------------------------------------------------#
    # Function tp shift indices to make a circular repressilator of 
    # arbitrary size - sign just add a negavtive sign for "repression"
    # Returns Matrix{Int64}
    #-------------------------------------------------------------------------#
    reg_rel = cyclic_permuations(regulators, sign = true)
    #-------------------------------------------------------------------------#
    # For the repressillator there is no need to have any other output
    # Except for TF binding since they bind and inhibit expression.
    # Build final grn Struct
    # Currently unmatable struct since once the GRN is set the rest happens
    # elsewhere. This just serves a behavior template
    #-------------------------------------------------------------------------#
    grn = GRN(regulatory_rel = reg_rel,
        regulators = regulators,
        strength_range = (0.99,1.0),
        regulator_strength = (-1) * strengths,
        messaging_output = nothing,
        metabolic_output = nothing,
        chromatin_remodelling = regulators,
        tf_binding = (-1) * regulators)
    return grn, gene_state
end

# function template_grn(regulator_strength::Vector{Float64};
#     n_regulators::Int64 = 1,
#     cellular_output::Int64 = 0,
#     metabolic_output::Int64 = 10,
#     chromatin_remodelling::Int64 = 0,
#     tf_binding::Int64 = 0,
#     overlap_range::Vector{Float64} = [0.0,0.0],
#     )::GRN
#     return 0
# end


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
    for g in eachindex(grns)
        grn,gene_state = repressilator(regulator_strength, gene_state)
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
    for g in eachindex(grns)
        grn,gene_state = repressilator(regulator_strength, gene_state)
        grn_set[grns[g]] = grn
    end
    sample.grn_set = grn_set
    sample.gene_state = gene_state
    return sample
end