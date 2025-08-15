using Plots
gr()  # or any other backend you like

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
function view_embeddings(sample::SampleState)
    out = transpose(Matrix(sample.out))
    pca = fit(PCA, out;maxoutdim=3)
end