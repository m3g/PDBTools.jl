import Plots: heatmap, scatter, cgrad
using SparseArrays: sparse, findnz

# the size of the plot should be proportional to the number 
# of residues in each dimension, with a 100 residues for 
# a 500px axis, by default
function _plot_size(map::ContactMap)
    nres1 = size(map.matrix, 1)
    nres2 = size(map.matrix, 2)
    xpixels = min(600, max(400, round(Int, 10 * nres1)))
    ypixels = min(600, max(400, round(Int, 10 * nres2)))
    return (xpixels, ypixels)
end

# Add space of not to colorbar title
_n(dmax) = dmax > 9.9 ? "\n" : ""

#
# The following "heatmap" function uses Plots.scatter under the hood, to deal 
# better with the non-attributed values of the sparse arrays.
#
"""
    heatmap(map::PDBTools.ContactMap; kwargs...)

Plot a contact map.

# Arguments

- `map::ContactMap`: the contact map to plot

All other arguments are default keywords of `Plots.heatmap` and can be adjusted to
customize the plot.

Most typical options to adjust are:

- `xstep`: the stride of the x-axis ticks 
- `ystep`: the stride of the y-axis ticks
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
ContactMap{Bool} of size (243, 12) with 17 contacts, threshold 4.0 and gap 0

julia> # plt = heatmap(map) # uncomment to plot

julia> map = contact_map(cA, cB; discrete=false) # distance map
ContactMap{Float32} of size (243, 12) with 17 contacts, threshold 4.0 and gap 0

julia> # plt = heatmap(map) # uncomment to plot
```

"""
function heatmap(
    map::ContactMap{T}; 
    plot_size=_plot_size(map),
    xstep=max(1, div(size(map.matrix, 1), 20)), 
    ystep=max(1, div(size(map.matrix, 2), 20)),
    xticks=PDBTools.residue_ticks(map.residues1; stride=xstep, serial=true),
    yticks=PDBTools.residue_ticks(map.residues2; stride=ystep, serial=true),
    xrotation=60,
    xlabel="residue",
    ylabel="residue",
    colorbar=ifelse(T <: Integer, :none, :right),
    colorbar_title=nothing,
    aspect_ratio=(last(plot_size)/first(plot_size))*(Base.size(map.matrix,1)/Base.size(map.matrix,2)),
    xlims=(0,size(map.matrix, 1)+1),
    ylims=(0,size(map.matrix, 2)+1),
    color=nothing,
    size=plot_size,
    framestyle=:box,
    grid=false,
    clims=nothing,
    margin=0.3Plots.Measures.cm,
    fontfamily="Computer Modern",
    markershape=:rect,
    markersize=2,
    markerstrokewidth=0,
    label="",
    kargs...
) where {T<:Real}
    i, j, invd = findnz(map.matrix)
    d = inv.(invd)
    ext = extrema(d)
    color = isnothing(color) ? 
        T == Bool ? cgrad([:white, :black], 2) : 
        (ext[1] < 0 ? :bwr : :grayC) : color
    clims = isnothing(clims) ? (1.1 * min(0, ext[1]), 1.1 * max(0, ext[2])) : clims
    colorbar_title = isnothing(colorbar_title) ?  
        ifelse(T <: Integer, nothing, "$(_n(ext[2]))distance (Å)") : 
        colorbar_tile
    return Plots.scatter(i, j; zcolor=d, colorbar, color,
        xlabel, ylabel, xticks, yticks, xrotation, label,
        colorbar_title, aspect_ratio, xlims, ylims,
        size, framestyle, grid, clims, margin, 
        markershape, markersize, markerstrokewidth,
        fontfamily, kargs...
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
    c_cont = contact_map(cA; discrete=false)
    plt = heatmap(c_cont)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
    c_cont.matrix .= rand.(Ref([-1,1])) .* c_cont.matrix 
    plt = heatmap(c_cont)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
end
