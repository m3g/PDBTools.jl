import Plots: plot

# Color of a text label drawn over a marker of color `c`: white for dark colors, black for
# light ones (relative luminance of the color).
_label_color(c) = 0.299 * Plots.red(c) + 0.587 * Plots.green(c) + 0.114 * Plots.blue(c) < 0.5 ? :white : :black

"""
    plot(wheel::HelicalWheel; kwargs...)

Plot a helical wheel, as computed by [`helical_wheel`](@ref), with each residue colored
by its value.

!!! note
    This function requires `using Plots`.

# Keyword arguments
- `normalize = true`: apply min-max normalization of the values to `[0, 1]`.
- `colormap = :tempo`: color map of the values.
- `colorbar_title = ""`: title of the color bar.
- `clims = nothing`: fixed limits of the color scale, e.g. `(0, 1)`.
  Use (with `normalize = false`) to compare several helices on the same scale.
- `show_moment = true`: draw the arrow of the hydrophobic moment.

The residues are placed at 100° from each other; after 18 residues (5 turns) a new ring
starts, with a larger radius. The wheel is rotated such that the hydrophobic moment points
down. The first and last residues are labeled with their residue numbers.

# Example

```julia-repl
julia> using PDBTools, Plots

julia> wheel = helical_wheel("LKKLLKLLKKLLKL", rand(14));

julia> plot(wheel; colorbar_title = "property")
```
"""
function plot(
    wheel::PDBTools.HelicalWheel;
    normalize::Bool = true,
    colormap = :tempo,
    colorbar_title::AbstractString = "",
    clims = nothing,
    show_moment::Bool = true,
)
    seq = collect(wheel.sequence)
    n = length(seq)
    μx, μy = wheel.μx, wheel.μy

    # Optional min-max normalization
    weights = copy(wheel.values)
    if normalize
        min_w, max_w = extrema(weights)
        weights = (weights .- min_w) ./ (max_w - min_w + 1e-9)
    end

    # Geometry: 100° per residue; new ring every 18 residues (5 turns)
    turn_residues = 18
    radius_step = 0.25
    angles = PDBTools._helical_wheel_angles(n)
    radii = [1 + div(i, turn_residues) * radius_step for i in 0:n-1]

    # Rotate the diagram so that the moment points down
    rotation_angle = atan(μy, μx) + π / 2
    rotated_angles = angles .- rotation_angle
    μx_rot = μx * cos(rotation_angle) + μy * sin(rotation_angle)
    μy_rot = -μx * sin(rotation_angle) + μy * cos(rotation_angle)

    x = cos.(rotated_angles) .* radii
    y = sin.(rotated_angles) .* radii

    # Limits of the color scale, used both for the markers and for the colors of their
    # labels (if all values are equal, the limits are widened, such that the markers
    # have the color of the center of the scale).
    lo, hi = isnothing(clims) ? extrema(weights) : clims
    if lo == hi
        lo, hi = lo - 0.5, hi + 0.5
    end

    gradient = Plots.cgrad(colormap)
    p = Plots.plot(x, y; lw = 1, lc = :gray10, alpha = 0.6, label = "", fontfamily = "Computer Modern")
    Plots.scatter!(p, x, y;
        marker_z = weights,
        markersize = 14,
        color = gradient,
        clims = (lo, hi),
        colorbar = :left,
        colorbar_title = colorbar_title,
        label = "",
    )

    # Residue labels: white on dark markers, black on light ones. The color of each
    # marker is that of its value in the color gradient, within the limits of the color scale.
    for i in 1:n
        t = clamp((weights[i] - lo) / (hi - lo), 0, 1)
        Plots.annotate!(p, x[i], y[i], Plots.text(seq[i], _label_color(gradient[t]), :center))
    end
    Plots.annotate!(p, x[1], y[1] - 0.16, Plots.text("$(first(wheel.residue_numbers))", :gray10, :center, 4))
    Plots.annotate!(p, x[end], y[end] - 0.16, Plots.text("$(last(wheel.residue_numbers))", :gray10, :center, 4))

    if show_moment
        Plots.quiver!(p, [0], [0]; quiver = ([μx_rot], [μy_rot]), color = :gray30, linewidth = 1, label = "")
    end

    lim = maximum(radii) + 0.5
    Plots.plot!(p;
        xlim = (-lim, lim), ylim = (-lim, lim),
        aspect_ratio = 1, grid = false, legend = false, frame = :box, ticks = [],
        dpi = 300,
    )

    return p
end

@testitem "plot(::HelicalWheel)" begin
    using PDBTools, LazyArtifacts
    import Plots
    ENV["GKSwstype"] = "nul"

    seq = "LKKLLKLLKKLLKL"
    n = length(seq)
    p = Plots.plot(helical_wheel(seq, collect(1.0:n)); colorbar_title = "CN")
    @test p isa Plots.Plot
    tmpfile = tempname() * ".png"
    Plots.savefig(p, tmpfile)
    @test isfile(tmpfile)
    # more than 18 residues (second ring), fixed color limits, no arrow
    p = Plots.plot(helical_wheel(repeat(seq, 2), rand(2n)); normalize = false, clims = (0, 1), show_moment = false)
    @test p isa Plots.Plot
    # all values equal (e.g. zero B-factors), with and without normalization
    @test Plots.plot(helical_wheel(seq, zeros(n))) isa Plots.Plot
    @test Plots.plot(helical_wheel(seq, zeros(n)); normalize = false) isa Plots.Plot
    # labels are white on dark markers and black on light ones
    label_color = Base.get_extension(PDBTools, :PlotsPlotting)._label_color
    gradient = Plots.cgrad(:tempo)
    @test label_color(gradient[0.0]) == :black
    @test label_color(gradient[1.0]) == :white
    @test label_color(Plots.RGB(0, 0, 0)) == :white
    @test label_color(Plots.RGB(1, 1, 1)) == :black

    # from a selection of atoms
    pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    atoms = select(pdb, at -> resnum(at) <= 18)
    p = Plots.plot(helical_wheel(atoms); colorbar_title = "B-factor")
    @test p isa Plots.Plot
end
