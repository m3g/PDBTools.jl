#
# Contact maps: heatmap(::ContactMap)
#

# Size of the figure proportional to the number of residues in each dimension, within
# 400 and 600 pt (plus the space of the color bar)
function _contact_map_size(map::ContactMap{T}) where {T}
    n1, n2 = Base.size(map.matrix)
    w = min(600, max(400, 10 * n1)) + (T <: Integer ? 0 : 80)
    h = min(600, max(400, 10 * n2))
    return (w, h)
end

"""
    heatmap(map::ContactMap; kwargs...)

Plot a contact map, as computed by [`contact_map`](@ref). Returns a `Makie.Figure`.

!!! note
    This function requires a Makie backend, e.g. `using CairoMakie`.

# Keyword arguments
- `colormap = nothing`: color map. If `nothing`: black for binary maps, `:grayC` for
  distance maps, and `:bwr` for maps with negative values (differences of maps: blue for
  contacts present only in the first map, red for those present only in the second).
- `clims = nothing`: limits of the color scale, as a tuple `(min, max)`.
- `colorbar`, `colorbar_title`: color bar, shown by default for distance maps, with the
  title `"distance (Å)"`. Not shown, by default, for binary maps and their differences.
- `xstep`, `ystep`: stride of the residue labels on the axes (default: about 20 labels per axis).
- `xlabel = "residue"`, `ylabel = "residue"`, `title = ""`: labels of the axes, and title.
- `size`, `font = nothing`, `fontsize = 14`: size of the figure (default: proportional to the
  number of residues), font (default: New Computer Modern, distributed with Makie) and font size.

# Example

```julia-repl
julia> using PDBTools, LazyArtifacts, CairoMakie

julia> ats = read_pdb(PDBTools.test_file("dimer.pdb"));

julia> map = contact_map(select(ats, "chain A"), select(ats, "chain B"); discrete = false);

julia> heatmap(map)
```
"""
function Makie.heatmap(map::ContactMap; size = _contact_map_size(map), font = nothing, fontsize = 14, kwargs...)
    fig = _figure(; size, font, fontsize)
    Makie.heatmap(fig[1, 1], map; kwargs...)
    return fig
end

function Makie.heatmap(gp::GridPos, map::ContactMap{T};
    colorbar::Bool = !(T <: Integer),
    colorbar_title::AbstractString = T <: Integer ? "" : "distance (Å)",
    xstep::Integer = max(1, div(Base.size(map.matrix, 1), 20)),
    ystep::Integer = max(1, div(Base.size(map.matrix, 2), 20)),
    xlabel::AbstractString = "residue",
    ylabel::AbstractString = "residue",
    title::AbstractString = "",
    kwargs...,
) where {T}
    n1, n2 = Base.size(map.matrix)
    gl = GridLayout(gp)
    ax = Axis(gl[1, 1];
        xlabel, ylabel, title,
        xticks = PDBTools.residue_ticks(map.residues1; stride = xstep, serial = true),
        yticks = PDBTools.residue_ticks(map.residues2; stride = ystep, serial = true),
        xticklabelrotation = π / 3,
        xgridvisible = false,
        ygridvisible = false,
    )
    plt = Makie.heatmap!(ax, map; kwargs...)
    limits!(ax, 0.5, n1 + 0.5, 0.5, n2 + 0.5)
    colorbar && Colorbar(gl[1, 2], plt; label = colorbar_title)
    return gl
end

function Makie.heatmap!(ax::Axis, map::ContactMap{T}; colormap = nothing, clims = nothing) where {T}
    # Dense matrix of distances, with NaN (transparent) where there is no contact. The
    # matrix stores the inverse of the distances (1 for the contacts of binary maps).
    n1, n2 = Base.size(map.matrix)
    m = fill(NaN, n1, n2)
    i, j, invd = findnz(map.matrix)
    d = isempty(invd) ? Float64[] : inv.(Float64.(invd))
    m[CartesianIndex.(i, j)] .= d
    ext = isempty(d) ? (0.0, 1.0) : extrema(d)
    if T == Bool
        colormap = isnothing(colormap) ? [:black, :black] : colormap
        clims = isnothing(clims) ? (0.0, 1.0) : clims
    else
        colormap = isnothing(colormap) ? (ext[1] < 0 ? :bwr : :grayC) : colormap
        clims = isnothing(clims) ? (1.1 * min(0, ext[1]), 1.1 * max(0, ext[2])) : clims
    end
    return Makie.heatmap!(ax, 1:n1, 1:n2, m; colormap, colorrange = clims, nan_color = :transparent)
end

@testitem "heatmap(::ContactMap)" begin
    using PDBTools, LazyArtifacts
    import CairoMakie
    tmpplot = tempname() * ".png"
    ats = read_pdb(PDBTools.test_file("dimer.pdb"))
    cA = select(ats, "chain A")
    cB = select(ats, "chain B")
    for map in (contact_map(cA, cB), contact_map(cA, cB; discrete = false), contact_map(cA))
        fig = CairoMakie.heatmap(map)
        @test fig isa CairoMakie.Figure
        CairoMakie.save(tmpplot, fig)
        @test isfile(tmpplot)
    end
    # difference of binary maps (ContactMap{Int}): no color bar by default
    cdiff = contact_map(cA; dmax = 8.0) - contact_map(cA)
    @test CairoMakie.heatmap(cdiff) isa CairoMakie.Figure
    # signed values (e.g. a difference of distance maps)
    c_cont = contact_map(cA; discrete = false)
    c_cont.matrix .= rand.(Ref([-1, 1])) .* c_cont.matrix
    @test CairoMakie.heatmap(c_cont; colorbar_title = "Δ distance (Å)") isa CairoMakie.Figure
    @test Tuple(CairoMakie.heatmap(c_cont; size = (300, 300)).scene.viewport[].widths) == (300, 300)
    # into a position of a figure, and into an axis
    fig = CairoMakie.Figure()
    @test CairoMakie.heatmap(fig[1, 1], contact_map(cA, cB)) isa CairoMakie.GridLayout
    @test CairoMakie.heatmap!(CairoMakie.Axis(fig[1, 2]), contact_map(cA, cB)) isa CairoMakie.AbstractPlot
end
