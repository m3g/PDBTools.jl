"""
    heatmap(map::HammerAitoffMap; kwargs...)

Plot a heat map of the values of a [`HammerAitoffMap`](@ref), as computed by
[`hammer_aitoff_map`](@ref).

!!! note
    This function requires `using Plots`.

# Keyword arguments
- `colormap=:tempo`: color map of the heat map.
- `clims=nothing`: color limits as a tuple `(min, max)`. If `nothing`, they are
  taken from the data.
- `contours=:none`: contour mode, one of
  - `:none`: no contours;
  - `:chain` (or `:chains`): draw one contour line around each chain;
  - `:group`: draw contours around groups of atoms, where the group of an
    atom is given by `groupfun(atom)` (for example, groups of chains, with
    `groupfun = at -> groups[chain(at)]`).
- `groupfun=chain`: function mapping an atom to its group (only used when
  `contours=:group`).
- `contour_color=:grey30`: color of the contour lines.
- `colorbar_title=""`, `title=""`: texts of the color bar and of the plot.

The plot can be saved with `savefig`.

# Example

```julia-repl
julia> using PDBTools, LazyArtifacts, Plots

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> map = hammer_aitoff_map(atoms, beta.(atoms));

julia> heatmap(map; contours = :chain, colorbar_title = "B-factor")
```
"""
function heatmap(
    map::PDBTools.HammerAitoffMap;
    colormap=:tempo,
    clims=nothing,
    contours::Symbol=:none,
    groupfun=PDBTools.chain,
    contour_color=:grey30,
    colorbar_title::AbstractString="",
    title::AbstractString="",
)
    if !(contours in (:none, :chain, :chains, :group))
        throw(ArgumentError("contours must be :none, :chain (or :chains) or :group, got $contours"))
    end

    p = Plots.contourf(
        map.xbin, map.ybin, map.values;
        color=colormap,
        clims=clims,
        xticks=[],
        yticks=[],
        lw=0.0,
        label="",
        colorbar_title=colorbar_title,
        title=title,
        framestyle=:box,
        fontfamily="Computer Modern",
    )

    if contours != :none
        # Label of each atom: its chain, or its group
        atoms = map.projection.atoms
        labels = contours == :group ? [string(groupfun(at)) for at in atoms] : [string(PDBTools.chain(at)) for at in atoms]

        # The color scale is shared by all series of the plot. To keep the
        # contour masks from changing it, the masks use only the lowest and
        # highest values of the color scale, and `clims` is passed to every
        # `contour!` call.
        lo, hi = isnothing(clims) ? PDBTools._finite_extrema(map.values) : clims

        # One contour per label: the mask is `hi` for atoms with that label
        # and `lo` for all the others, and the contour is drawn at the border.
        # The bins of the masks are the same of the map, because they are
        # computed from the same projected points.
        nbins = length(map.xbin)
        for label in unique(labels)
            mask_values = [l == label ? hi : lo for l in labels]
            xb, yb, mask = PDBTools._binned_statistic_2d(map.projection.x, map.projection.y, mask_values; nbins, statistic=:mean)
            # Bins without atoms (NaN) do not belong to the group: contour lines are
            # not drawn correctly across NaN values
            mask[isnan.(mask)] .= lo
            Plots.contour!(
                p, xb, yb, mask;
                levels=1,
                clims=(lo, hi),
                color=contour_color,
                lw=0.5,
                label="",
            )
        end
    end

    return p
end

@testitem "heatmap(::HammerAitoffMap)" begin
    using PDBTools, LazyArtifacts
    import Plots
    ENV["GKSwstype"] = "nul"
    atoms = read_mmcif(PDBTools.test_file("6CO8.cif"), "protein")
    n = length(atoms)
    map = hammer_aitoff_map(atoms, beta.(atoms); nbins = 30)
    for contours in (:none, :chain)
        tmpfile = tempname() * ".png"
        p = Plots.heatmap(map; contours)
        @test p isa Plots.Plot
        Plots.savefig(p, tmpfile)
        @test isfile(tmpfile)
    end
    map = hammer_aitoff_map(atoms, ones(n); nbins = 30)
    p = Plots.heatmap(map; contours = :group, groupfun = at -> chain(at) == "A" ? "A" : "other")
    @test p isa Plots.Plot
    @test_throws ArgumentError Plots.heatmap(map; contours = :residue)
end
