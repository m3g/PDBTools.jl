#
# Hammer-Aitoff maps: heatmap(::HammerAitoffMap)
#

"""
    heatmap(map::HammerAitoffMap; kwargs...)

Plot a heat map of the values of a [`HammerAitoffMap`](@ref), as computed by
[`hammer_aitoff_map`](@ref). Returns a `Makie.Figure`.

!!! note
    This function requires a Makie backend, e.g. `using CairoMakie`.

# Keyword arguments
- `colormap = :tempo`: color map of the heat map.
- `clims = nothing`: color limits as a tuple `(min, max)`. If `nothing`, they are
  taken from the data.
- `contours = :none`: contour mode, one of
  - `:none`: no contours;
  - `:chain` (or `:chains`): draw one contour line around each chain;
  - `:group`: draw contours around groups of atoms, where the group of an
    atom is given by `groupfun(atom)` (for example, groups of chains, with
    `groupfun = at -> groups[chain(at)]`).
- `groupfun = chain`: function mapping an atom to its group (only used when
  `contours = :group`).
- `contour_color = :grey30`: color of the contour lines.
- `colorbar = true`, `colorbar_title = ""`, `title = ""`: color bar and title of the plot.
- `size = (700, 400)`, `font = nothing`, `fontsize = 14`: size of the figure, font
  (default: New Computer Modern, distributed with Makie) and font size.

# Example

```julia-repl
julia> using PDBTools, LazyArtifacts, CairoMakie

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> map = hammer_aitoff_map(atoms, beta.(atoms));

julia> heatmap(map; contours = :chain, colorbar_title = "B-factor")
```
"""
function Makie.heatmap(map::HammerAitoffMap; size = (700, 400), font = nothing, fontsize = 14, kwargs...)
    fig = _figure(; size, font, fontsize)
    Makie.heatmap(fig[1, 1], map; kwargs...)
    return fig
end

function Makie.heatmap(gp::GridPos, map::HammerAitoffMap;
    colorbar::Bool = true,
    colorbar_title::AbstractString = "",
    title::AbstractString = "",
    kwargs...,
)
    gl = GridLayout(gp)
    ax = Axis(gl[1, 1]; title, aspect = DataAspect())
    hidedecorations!(ax)
    plt = Makie.heatmap!(ax, map; kwargs...)
    colorbar && Colorbar(gl[1, 2], plt; label = colorbar_title)
    return gl
end

function Makie.heatmap!(ax::Axis, map::HammerAitoffMap;
    colormap = :tempo,
    clims = nothing,
    contours::Symbol = :none,
    groupfun = PDBTools.chain,
    contour_color = :grey30,
)
    if !(contours in (:none, :chain, :chains, :group))
        throw(ArgumentError("contours must be :none, :chain (or :chains) or :group, got $contours"))
    end
    lo, hi = isnothing(clims) ? PDBTools._finite_extrema(map.values) : clims
    if lo == hi # all values equal
        lo, hi = lo - 0.5, hi + 0.5
    end
    # values[j, i] is the value of the bin with center (xbin[i], ybin[j])
    plt = contourf!(ax, map.xbin, map.ybin, permutedims(map.values);
        levels = range(lo, hi; length = 16),
        colormap,
        extendlow = :auto,
        extendhigh = :auto,
    )

    if contours != :none
        # Label of each atom: its chain, or its group
        atoms = map.projection.atoms
        labels = contours == :group ? [string(groupfun(at)) for at in atoms] : [string(PDBTools.chain(at)) for at in atoms]
        # One contour per label: the mask is 1 for atoms with that label and 0 for all
        # the others, and the contour is drawn at 0.5. The bins of the masks are the same
        # of the map, because they are computed from the same projected points.
        nbins = length(map.xbin)
        for label in unique(labels)
            mask_values = [l == label ? 1.0 : 0.0 for l in labels]
            xb, yb, mask = PDBTools._binned_statistic_2d(map.projection.x, map.projection.y, mask_values; nbins, statistic = :mean)
            # Bins without atoms (NaN) do not belong to the group
            mask[isnan.(mask)] .= 0.0
            contour!(ax, xb, yb, permutedims(mask); levels = [0.5], color = contour_color, linewidth = 0.75)
        end
    end
    return plt
end

@testitem "heatmap(::HammerAitoffMap)" begin
    using PDBTools, LazyArtifacts
    import CairoMakie
    atoms = read_mmcif(PDBTools.test_file("6CO8.cif"), "protein")
    n = length(atoms)
    map = hammer_aitoff_map(atoms, beta.(atoms); nbins = 30)
    for contours in (:none, :chain)
        fig = CairoMakie.heatmap(map; contours, colorbar_title = "B-factor")
        @test fig isa CairoMakie.Figure
        tmpfile = tempname() * ".png"
        CairoMakie.save(tmpfile, fig)
        @test isfile(tmpfile)
    end
    # constant values; contours around groups; no color bar
    map1 = hammer_aitoff_map(atoms, ones(n); nbins = 30)
    fig = CairoMakie.heatmap(map1; contours = :group, groupfun = at -> chain(at) == "A" ? "A" : "other", colorbar = false)
    @test fig isa CairoMakie.Figure
    @test Tuple(CairoMakie.heatmap(map; size = (500, 300)).scene.viewport[].widths) == (500, 300)
    @test_throws ArgumentError CairoMakie.heatmap(map; contours = :residue)
    @test_throws ArgumentError CairoMakie.heatmap(map; size = (0, 300))
    # into a position of a figure, and into an axis
    fig = CairoMakie.Figure()
    @test CairoMakie.heatmap(fig[1, 1], map) isa CairoMakie.GridLayout
    @test CairoMakie.heatmap!(CairoMakie.Axis(fig[2, 1]), map; clims = (0, 50)) isa CairoMakie.AbstractPlot
end
