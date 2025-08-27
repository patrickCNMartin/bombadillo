using Plots
function plot_tissue_2d_old(tissue::Tissue)
    coords = tissue.coordinates
    labels = tissue.cell_labels

    if coords === nothing || labels === nothing
        error("Tissue must have both `coordinates` and `cells` to plot.")
    end

    
    x = getindex.(coords,1)
    y = getindex.(coords,2)

    Plots.scatter(
        x,
        y, 
        zcolor = labels  
    )
end


#-----------------------------------------------------------------------------#
# New Viz function for cell states 
#-----------------------------------------------------------------------------#
function plot_tissue(sample::SampleState,
    labels::Symbol = :cell_types;
    dim = 2)
    gr() # to be updated when I want to save plots latter
    tissue = sample.tissue
    coords = tissue.coordinates
    labs = getfield(tissue, labels)
    if isnothing(labs)
        throw(DomainError("Requested label is not present in sample"))
    end
    x = getindex.(coords,1)
    y = getindex.(coords,2)

    if dim == 2
        Plots.scatter(x,y,group = labs)
    elseif dim == 3
        z = getindex.(coords,3)
        Plots.scatter(x,y,z,group = labs)
    else
        throw(ArgumentError("dim should either be 2 or 3"))
    end
end


#-----------------------------------------------------------------------------#
# Viz functions for spatial data stuff
# these are all very rough functions that I might not keep in the long run
# I just want to see what it looks like at test a few ideas 
#-----------------------------------------------------------------------------#
using Plots
using Statistics

function view_gene_cycle(
    sample::SampleState,
    gene_pull::Vector{DataFrame};
    sort::Bool = true,
    merge_cells::Bool = true,
    norm_rank::Bool = true
)
    gr()
    panels = ceil(Int, sqrt(length(gene_pull)))
    # making copy of gene pull to avoid any change at global level
    gene_shift = deepcopy(gene_pull)
    if merge_cells
        gene_vec = Vector{Vector{Float64}}(undef, length(gene_shift))
        min_y = 2000 # Use a starting value with maximum potential rnak shift
        max_y = 0.0
        n_cells = size(gene_shift[1],1)
        n_cycles = size(gene_shift[1],2)
        for (i, g) in enumerate(gene_shift)
            d = zeros(Float64, n_cells, n_cycles)
            for t in 2:n_cycles
                @views d[:, t] .= g[:, t] .- g[:, t-1]
            end
            diff_mat = g .- d
            nz = [sum(x) for x in eachrow(diff_mat)]
            nz = findall(x -> x !=0, nz)
            mean_shift = [mean(x) for x in eachcol(diff_mat)]
            min_y = minimum([min_y, minimum(mean_shift[mean_shift .!= 0])])
            max_y = maximum([max_y, maximum(mean_shift[mean_shift .!= 0])])
            gene_vec[i] = mean_shift
            
        end

        P = plot(ylims=(min_y, max_y), xlabel="Cycle", ylabel="Mean Rank Change", title="Gene Rank Shifts")
        for (p, g) in enumerate(gene_vec)
            x = 1:length(g)
            plot!(P, x, g, label="Gene $p")
        end
        display(P)

    else
        if sort
            for g in gene_shift
                numeric_cols = names(g, Number)
                row_vars = [var(collect(skipmissing(r))) for r in eachrow(g[:, numeric_cols])]
                g.row_variance = row_vars
                sort!(g, :row_variance)
                select!(g, Not(:row_variance))
            end
        end

        plot_array = Vector{Any}()
        for (p, g) in enumerate(gene_shift)
            numeric_cols = names(g, Number)
            local_plot = heatmap(Matrix(g[:, numeric_cols]),
                                 title="Gene $p", colorbar=true)
            push!(plot_array, local_plot)
        end

        plot(plot_array..., layout=(panels, panels))
    end
end


function view_tissue_cycle(
    sample::SampleState,
    gene_pull::Vector{DataFrame},
    output_path::String;
    dim::Int = 2,
    fps::Int = 10,
    markersize::Real = 3)
    plotly()
    try
        coordinates = sample.tissue.coordinates
        x = getindex.(coordinates, 1)
        y = getindex.(coordinates, 2)
        z = getindex.(coordinates, 3)

        dim ∈ (2, 3) || throw(ArgumentError("dim should be either 2 or 3"))

        num_time = size(gene_pull[1], 2)
        for p in gene_pull[2:end]
            size(p, 2) == num_time || throw(ArgumentError("All DataFrames must have the same number of columns (time points)"))
        end

        n_genes = length(gene_pull)
        n_cols = ceil(Int, sqrt(n_genes))
        n_rows = ceil(Int, n_genes / n_cols)
        layout = (n_rows, n_cols)

        n_cells = nrow(gene_pull[1])
        diff_mats = Vector{Matrix{Float64}}(undef, n_genes)
        global_max = 0.0
        for g in 1:n_genes
            diff_mat = delta_shift(gene_pull[g])
            diff_mats[g] = diff_mat
            m = maximum(abs, diff_mat)
            if isfinite(m) && m > global_max
                global_max = m
            end
        end
        global_max = global_max == 0 ? 1.0 : global_max
        clims = (-global_max, global_max)
        diverge = cgrad(["#2f4261", "#f5f3e4", "#9c2222"])

        anim = @animate for t in 1:num_time
            subplots = Plots.Plot[]
            for g in 1:n_genes
                vals = @views diff_mats[g][:, t]
                alpha_vals = 0.25 .+ 0.75 .* (abs.(vals) ./ global_max)   # scale 0 → 1
                plt = if dim == 2
                    scatter(x, y;
                        zcolor = vals,
                        c = diverge,
                        clims = clims,
                        alpha = alpha_vals,
                        legend = false,
                        colorbar = true,
                        ms = markersize,
                        markerstrokewidth = 0,
                        title = "Gene $(g) Δ")
                else
                    scatter(x, y, z;
                        zcolor = vals,
                        c = diverge,
                        clims = clims,
                        alpha = alpha_vals,
                        legend = false,
                        colorbar = true,
                        ms = markersize,
                        markerstrokewidth = 0,
                        title = "Gene $(g) Δ")
                end
                push!(subplots, plt)
            end
            plot(subplots...; layout = layout, size = (420*n_cols, 420*n_rows))
        end

        gif(anim, output_path; fps = fps)
    finally
        gr()
    end
end


function delta_shift(pull::DataFrame)::Matrix
    n_cells = size(pull, 1)
    num_time = size(pull, 2)
    delta_mat = zeros(Float64, n_cells, num_time)
    for t in 2:num_time
        @views delta_mat[:, t] .= pull[:, t] .- pull[:, t-1]
    end
    return delta_mat
end