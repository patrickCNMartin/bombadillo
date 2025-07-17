function generate_grn(
    n_grns:: Int64 = 100,
    tamplate:: String = "random")

end



function repressilator(
    regulator_strength::Vector{Float64};
    n_regulators::Int64 = 3)::GRN
    regulators, strengths = compute_grn_overlaps(
        regulator_strength,
        overlap_range = [0.0,0.0],
        stength_range = [-1.0,-0.8],
        g = n_regulators)
    
    


    cellular_output = zeros(n_regulators)
    metabolic_output = zeros(n_regulators)
    chromatin_remodelling = zeros(n_regulators)
    tf_binding = regulators
    grn = GRN(regulatory_rel = regulatory_rel,
        regulators = regulators,
        regulator_strength = strengths,
        cellular_output = cellular_output,
        metabolic_output = metabolic_output,
        chromatin_remodelling = chromatin_remodelling,
        tf_binding = tf_binding)
    return grn
end


function compute_grn_overlaps(
    regulator_strength::Vector{Float64};
    overlap_range::Vector{Float64} = [0.0,0.0],
    stength_range::Vector{Float64} = [0.5,1.0],
    g::Int64 = 3)::Tuple{Vector{Int64},Vector{Float64}}
    if overlap_range[1] == 0.0 && overlap_range[2] == 0.0
        locs = findall(regulator_strength .== 0.0)
        regulators = rand(locs, g)
        strenghts = rand(Uniform(stength_range[1]:stength_range[2]),g) 
    else
        over_g = rand(Int(g * overlap_range[1]):Int(g * overlap_range[2]))
        new_g = g - over_g
        over_locs = findall(regulator_strength .> 0)
        new_locs = findall(regulator_strength .== 0)
        regulators = [rand(over_locs,over_g);rand(new_locs,new_g)]
        strenghts = rand(Uniform(stength_range[1]:stength_range[2]), over_g + new_g)
    end
    return regulators, strenghts
end