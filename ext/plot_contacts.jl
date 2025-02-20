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

"""
    heatmap(map::ContactMap; kwargs...)

Plot a contact map.

# Arguments

- `map::ContactMap`: the contact map to plot

All other arguments are default keywords of `Plots.heatmap` and can be adjusted to
customize the plot.

Most typical options to adjust are:

- `color`: the color palette to use (default: `:grayC` for distances, `:Greys_9` for binary maps)
- `clims`: the range of the color scale.
- `colorbar_title`: the title of the colorbar. Default: "distance (Å)" for distances, no title for binary maps.

# Example 

```jldoctest
julia> using PDBTools, Plots

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> cA = select(ats, "chain A");

julia> cB = select(ats, "chain B");

julia> map = contact_map(cA, cB)
ContactMap{Bool} of size (243, 12), with threshold 4.0 and gap 0

julia> plt = heatmap(map) # produced the figure

julia> map = contact_map(cA, cB; discrete=false) # distance map
ContactMap{Float32} of size (243, 12), with threshold 4.0 and gap 0

julia> plt = heatmap(map) # produces the figure
"""
function heatmap(::ContactMap) end

# heatmap for distance (quantitative) maps
function heatmap(
    map::ContactMap{Union{Missing,T}}; 
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
    clims=(0,1.1) .* extrema(skipmissing(map.matrix)),
    margin=0.5Plots.Measures.cm,
    kargs...
) where {T<:Real}
    return heatmap(transpose(map.matrix); 
        xlabel, ylabel, xticks, yticks, xrotation,
        colorbar_title, color, aspect_ratio, xlims, ylims,
        size, framestyle, grid, clims, margin,
        kargs...
    )
end

# heatmap for binary (discrete) maps
function heatmap(
    map::ContactMap{Union{Missing, Bool}}; 
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

@testitem "contact plots" begin
    using PDBTools
    using Plots
    tmpplot = tempname()*".png"
    ats = read_pdb(PDBTools.DIMERPDB)
    cA = select(ats, "chain A")
    cB = select(ats, "chain B")
    plt = heatmap(contact_map(cA, cB))
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
    plt = heatmap(contact_map(cA, cB; discrete=false))
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
    plt = heatmap(contact_map(cA))
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
    plt = heatmap(contact_map(cA; discrete=false))
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
end
