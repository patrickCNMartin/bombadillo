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

function add_genes(
    tissue::Tissue,
    cells::BaseCells,
    log_expression_range::NamedTuple = (mean = 2, variance = 1),
    dispersion::Float64 = 1.0,
    scale_factor::Float64 = 2.0,
    zi_ratio::Float64 = 0.2)
    n_genes = cells.n_genes
    n_cells = length(tissue.cell_labels)
    cells = last.(tissue.ecosystems)
    count_matrix = Matrix{Float64}(undef,n_genes, n_cells)

    
    gene_sample = generate_gene_sample(
        log_expression_range,
        n_genes,
        dispersion,
        scale_factor,
        zi_ratio)
    
    
    for cell in eachindex(cells)
        local_noise = rand.(Uniform(0.975, 1.025), n_genes)
        noisy_profile = round.(Int, gene_sample .* local_noise)
        count_matrix[:,cell] = noisy_profile[Int.(cells[cell])]
    end
    return count_matrix
    
end

#-----------------------------------------------------------------------------#
# New gene state - Pute simulated for now.
# we will change it later to create "synthetic data"
#-----------------------------------------------------------------------------#
function initialize_genes(
    n_genes::Int64 = 2000;
    max_leak::Int64 = 10,
    max_decay::Int64 = 100,
    translation_efficiency::Vector{Float64} = [0.5, 1.0])::GeneState
    #-------------------------------------------------------------------------#
    # build initial gene set
    #-------------------------------------------------------------------------#
    genes = string.("gene_",1:n_genes)
    regulator_strength = zeros(Float64, n_genes)
    saturation_rank = rand([n_genes / 2, n_genes],n_genes)
    leak_rate = rand([0,max_leak],n_genes)
    decay_rate = rand([max_leak, max_decay],n_genes) .* (-1)
    translation_efficiency = rand(
        Uniform(translation_efficiency[1],translation_efficiency[2]), n_genes)
    gene_state = GeneState(n_genes = n_genes,
        genes = genes,
        saturation_rank = saturation_rank,
        leak_rate = leak_rate,
        decay_rate = decay_rate,
        translation_efficiency = translation_efficiency,
        regulator_strength = regulator_strength)
   
    return gene_state
end