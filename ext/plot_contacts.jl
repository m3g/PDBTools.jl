import Plots: heatmap, cgrad

# the size of the plot should be proportional to the number 
# of residues in each dimension, with a 100 residues for 
# a 500px axis, by default
function _plot_size(map::ContactMap)
    nres1 = size(map.matrix, 1)
    nres2 = size(map.matrix, 2)
    xpixels = min(600, max(300, round(Int,500 * nres1 / 100)))
    ypixels = min(600, max(300, round(Int,500 * nres2 / 100)))
    return (xpixels, ypixels)
end

function heatmap(
    map::ContactMap{<:Real}; 
    plot_size=_plot_size(map),
    xstep=max(1, div(size(map.matrix, 1), 20)), 
    ystep=max(1, div(size(map.matrix, 2), 20)),
    xticks=PDBTools.residue_ticks(map.residues1; stride=xstep, serial=true),
    yticks=PDBTools.residue_ticks(map.residues2; stride=ystep, serial=true),
    xrotation=60,
    xlabel="residue",
    ylabel="residue",
    colorbar_title="distance (Å)",
    aspect_ratio=(last(plot_size)/first(plot_size))*(Base.size(map.matrix,1)/Base.size(map.matrix,2)),
    xlims=(1,size(map.matrix, 1)),
    ylims=(1,size(map.matrix, 2)),
    color=:grayC,
    size=plot_size,
    framestyle=:box,
    grid=false,
    clims=(1,1.5) .* extrema(skipmissing(map.matrix)),
    margin=0.5Plots.Measures.cm,
    kargs...
)
    return heatmap(transpose(map.matrix); 
        xlabel, ylabel, xticks, yticks, xrotation,
        colorbar_title, color, aspect_ratio, xlims, ylims,
        size, framestyle, grid, clims, margin,
        kargs...
    )
end

function heatmap(
    map::ContactMap{<:Bool}; 
    plot_size=_plot_size(map),
    xstep=max(1, div(size(map.matrix, 1), 20)), 
    ystep=max(1, div(size(map.matrix, 2), 20)),
    xticks=PDBTools.residue_ticks(map.residues1; stride=xstep, serial=true),
    yticks=PDBTools.residue_ticks(map.residues2; stride=ystep, serial=true),
    xrotation=60,
    xlabel="residue",
    ylabel="residue",
    xlims=(1,size(map.matrix, 1)),
    ylims=(1,size(map.matrix, 2)),
    color=:Greys_9,
    size=plot_size,
    aspect_ratio=(last(plot_size)/first(plot_size))*(Base.size(map.matrix,1)/Base.size(map.matrix,2)),
    framestyle=:box,
    grid=false,
    clims=(0,1),
    colorbar=:none,
    margin=0.5Plots.Measures.cm,
    kargs...
)
    return heatmap(transpose(map.matrix); 
        xlabel, ylabel, xticks, yticks, xrotation,
        color, colorbar, xlims, ylims, aspect_ratio,
        size, framestyle, grid, clims, margin,
        kargs...
    )
end
