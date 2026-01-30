import Plots: heatmap, cgrad
using SparseArrays: sparse, nonzeros, findnz

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

# If the number of contacts is too large, abort and warn the user for the 
# use of the `unsafe` option
function _unsafe(unsafe, map; n_unsafe=10^8)
    n = length(map.matrix)
    if !unsafe & (n > n_unsafe)
        throw(ArgumentError("""\n
            The size of the contact matrix ($(size(map.matrix))) is too large.
            Plotting it may explode the memory of your computer.

            If you are sure you want to plot it, use the `unsafe` option:
            
                heatmap(map; unsafe=true)
            
            at your own risk.

        """))
    end
end

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
ContactMap{Bool} of size (243, 12), with threshold 4.0 and gap 0

julia> # plt = heatmap(map) # uncomment to plot

julia> map = contact_map(cA, cB; discrete=false) # distance map
ContactMap{Float32} of size (243, 12), with threshold 4.0 and gap 0

julia> # plt = heatmap(map) # uncomment to plot
```

"""
function heatmap(::ContactMap) end

# heatmap for distance (quantitative) maps
function heatmap(
    map::ContactMap{T}; 
    unsafe=false, n_unsafe=10^8,
    plot_size=_plot_size(map),
    xstep=max(1, div(size(map.matrix, 1), 20)), 
    ystep=max(1, div(size(map.matrix, 2), 20)),
    xticks=PDBTools.residue_ticks(map.residues1; stride=xstep, serial=true),
    yticks=PDBTools.residue_ticks(map.residues2; stride=ystep, serial=true),
    xrotation=60,
    xlabel="residue",
    ylabel="residue",
    colorbar=ifelse(T <: Integer, :none, :right),
    colorbar_title=ifelse(T <: Integer, nothing, "\ndistance (Å)"),
    aspect_ratio=(last(plot_size)/first(plot_size))*(Base.size(map.matrix,1)/Base.size(map.matrix,2)),
    xlims=(1,size(map.matrix, 1)),
    ylims=(1,size(map.matrix, 2)),
    color=nothing,
    size=plot_size,
    framestyle=:box,
    grid=false,
    clims=nothing,
    margin=0.3Plots.Measures.cm,
    fontfamily="Computer Modern",
    kargs...
) where {T<:Real}
    _unsafe(unsafe, map; n_unsafe)
    i, j, invd = findnz(map.matrix)
    d = inv.(invd)
    ext = extrema(d)
    distance_matrix = Matrix{Union{Missing, T}}(undef, Base.size(map.matrix))
    distance_matrix .= missing
    for iel in eachindex(i, j, d)
        distance_matrix[i[iel],j[iel]] = d[iel]
    end
    color = isnothing(color) ? (ext[1] < 0 ? :bwr : :grayC) : color
    clims = isnothing(clims) ? (1.1 * min(0, ext[1]), 1.1 * max(0, ext[2])) : clims
    return heatmap(transpose(distance_matrix); 
        xlabel, ylabel, xticks, yticks, xrotation,
        colorbar_title, color, aspect_ratio, xlims, ylims,
        size, framestyle, grid, clims, margin, colorbar, 
        fontfamily, kargs...
    )
end

# heatmap for binary (discrete) maps
function heatmap(
    map::ContactMap{Bool}; 
    unsafe=false, n_unsafe=10^8,
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
    margin=0.3Plots.Measures.cm,
    fontfamily="Computer Modern",
    kargs...
)
    _unsafe(unsafe, map; n_unsafe)
    return heatmap(transpose(map.matrix); 
        xlabel, ylabel, xticks, yticks, xrotation,
        color, colorbar, xlims, ylims, aspect_ratio,
        size, framestyle, grid, clims, margin, 
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

    # Test unsafe error
    cm = contact_map(cA)
    @test_throws "risk" heatmap(cm; n_unsafe=100)
    plt = heatmap(cm; unsafe=true, n_unsafe=100)
    savefig(plt, tmpplot)
    rm(tmpplot; force=true)
end
