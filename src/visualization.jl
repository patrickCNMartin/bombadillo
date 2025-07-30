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
function plot_tissue(tissue::TissueState,
    labels::Symbol = :cell_types;
    dim = 2)
    coords = tissue.coordinates
    labs = getfield(tissue, labels)
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