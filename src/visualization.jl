using Plots
gr()  # or any other backend you like

function plot_tissue_2d(tissue::Tissue)
    coords = tissue.coordinates
    labels = collect(values(tissue.cell_labels))

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
