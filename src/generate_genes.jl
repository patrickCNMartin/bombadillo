#-----------------------------------------------------------------------------#
# New gene state - Pute simulated for now.
# we will change it later to create "synthetic data"
#-----------------------------------------------------------------------------#
using StatsBase
function initialize_genes(
    n_genes::Int64 = 2000;
    leak_range::Tuple{Float64,Float64} = (0.0,20.0),
    decay_range::Tuple{Float64,Float64} = (0.5,0.8),
    saturation::Float64 = 0.1,
    translation_efficiency::Vector{Float64} = [0.8, 1.0])::GeneState
    #-------------------------------------------------------------------------#
    # build initial gene set
    #-------------------------------------------------------------------------#
    genes = string.("gene_",1:n_genes)
    regulator_strength = zeros(Float64, n_genes)
    remodeler_strength = zeros(Float64, n_genes)
    leak_rate = rand(Uniform(leak_range[1],leak_range[2]),n_genes)
    decay_rate = rand(Uniform(decay_range[1],decay_range[2]),n_genes)
    saturation = rand(Uniform(1,n_genes * saturation),n_genes)
    translation_efficiency = rand(
        Uniform(translation_efficiency[1],translation_efficiency[2]), n_genes)
    gene_state = GeneState(n_genes = n_genes,
        genes = genes,
        saturation = saturation,
        leak_rate = leak_rate,
        decay_rate = decay_rate,
        translation_efficiency = translation_efficiency,
        regulator_strength = regulator_strength,
        remodeler_strength = remodeler_strength)
   
    return gene_state
end


using Distributions
function generate_gene_sample(
    log_expression_range::NamedTuple,
    n_genes::Int64 = 2000,
    dispersion::Float64= 1.0,
    scale_factor::Float64 = 2.0,
    zi_ratio::Float64 = 0.2)
    μ = log_expression_range.mean
    σ = log_expression_range.variance

    log_dist = LogNormal(μ, σ)
    μ_genes = clamp.(rand(log_dist, n_genes),1.0, Inf)
    μ_genes = μ_genes .* scale_factor
    counts = Vector{Int64}(undef, n_genes)
    for g in eachindex(μ_genes)
        μ = μ_genes[g]
        r = 1 / dispersion
        p = r / (r + μ)
        nb = NegativeBinomial(r, p)
        exp = rand(nb)
        counts[g] = exp
    end
    if zi_ratio > 0.0
        dropped = round(Int,length(counts) * zi_ratio)
        zi = rand(1:n_genes, dropped)
        counts[zi] .= 0
    end

    return sort(counts)
end


function add_counts(
    sample::SampleState,
    log_expression_range::NamedTuple = (mean = 2, variance = 1),
    dispersion::Float64 = 1.0,
    scale_factor::Float64 = 2.0,
    zi_ratio::Float64 = 0.2,
    layer::Symbol = :rna_state)::SampleState
    #-------------------------------------------------------------------------#
    # Pull everything out
    #-------------------------------------------------------------------------#
    n_cells = sample.n_cells
    n_genes = sample.n_genes
    cells = sample.cells
    count_matrix = Matrix{Float64}(undef,n_genes, n_cells)

    #-------------------------------------------------------------------------#
    # Create one distribution to be sampled
    #-------------------------------------------------------------------------#
    gene_sample = generate_gene_sample(
        log_expression_range,
        n_genes,
        dispersion,
        scale_factor,
        zi_ratio)
    
    #-------------------------------------------------------------------------#
    # Now we sample and put into the matrix
    #-------------------------------------------------------------------------#
    for cell in eachindex(cells)
        state = ranks_from_p(getfield(cells[cell],layer))
        noisy_profile = reverse(gene_sample)
        count_matrix[:,cell] = noisy_profile[state]
    end
   
    sample.out = count_matrix
    return sample
end


