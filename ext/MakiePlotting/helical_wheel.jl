#
# Helical wheels: plot(::HelicalWheel)
#

# Color of a text label drawn over a marker of color `c`: white for dark colors, black for
# light ones (relative luminance of the color).
function _label_color(c)
    c = Makie.to_color(c)
    return 0.299 * c.r + 0.587 * c.g + 0.114 * c.b < 0.5 ? :white : :black
end

"""
    plot(wheel::HelicalWheel; kwargs...)

Plot a helical wheel, as computed by [`helical_wheel`](@ref), with each residue colored
by its value. Returns a `Makie.Figure`.

!!! note
    This function requires a Makie backend, e.g. `using CairoMakie`.

# Keyword arguments
- `normalize = true`: apply min-max normalization of the values to `[0, 1]`.
- `colormap = :tempo`: color map of the values.
- `clims = nothing`: fixed limits of the color scale, e.g. `(0, 1)`.
  Use (with `normalize = false`) to compare several helices on the same scale.
- `show_moment = true`: draw the arrow of the hydrophobic moment.
- `colorbar = true`, `colorbar_title = ""`, `title = ""`: color bar and title of the plot.
- `size = (500, 420)`, `font = nothing`, `fontsize = 14`: size of the figure, font
  (default: New Computer Modern, distributed with Makie) and font size.

The residues are placed at 100° from each other; after 18 residues (5 turns) a new ring
starts, with a larger radius. The wheel is rotated such that the hydrophobic moment points
down. The first and last residues are labeled with their residue numbers.

# Example

```julia-repl
julia> using PDBTools, CairoMakie

julia> wheel = helical_wheel("LKKLLKLLKKLLKL", rand(14));

julia> plot(wheel; colorbar_title = "property")
```
"""
function Makie.plot(wheel::HelicalWheel; size = (500, 420), font = nothing, fontsize = 14, kwargs...)
    fig = _figure(; size, font, fontsize)
    Makie.plot(fig[1, 1], wheel; kwargs...)
    return fig
end

function Makie.plot(gp::GridPos, wheel::HelicalWheel;
    colorbar::Bool = true,
    colorbar_title::AbstractString = "",
    title::AbstractString = "",
    kwargs...,
)
    gl = GridLayout(gp)
    ax = Axis(gl[1, 1]; title, aspect = DataAspect())
    hidedecorations!(ax)
    plt = Makie.plot!(ax, wheel; kwargs...)
    colorbar && Colorbar(gl[1, 2], plt; label = colorbar_title)
    return gl
end

function Makie.plot!(ax::Axis, wheel::HelicalWheel;
    normalize::Bool = true,
    colormap = :tempo,
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

    gradient = cgrad(colormap)
    lines!(ax, x, y; color = (:gray10, 0.6), linewidth = 1)
    plt = scatter!(ax, x, y;
        color = weights,
        colormap = gradient,
        colorrange = (lo, hi),
        markersize = 28,
        strokewidth = 0.5,
        strokecolor = :gray30,
    )
    # Residue labels: white on dark markers, black on light ones. The color of each
    # marker is that of its value in the color gradient, within the limits of the color scale.
    label_colors = [_label_color(gradient[clamp((w - lo) / (hi - lo), 0, 1)]) for w in weights]
    text!(ax, x, y; text = string.(seq), color = label_colors, align = (:center, :center), fontsize = 14)
    text!(ax, [x[1], x[end]], [y[1], y[end]] .- 0.17;
        text = string.([first(wheel.residue_numbers), last(wheel.residue_numbers)]),
        color = :gray10, align = (:center, :center), fontsize = 8)

    if show_moment
        arrows2d!(ax, [Point2f(0, 0)], [Vec2f(μx_rot, μy_rot)]; color = :gray30, shaftwidth = 1.5)
    end

    lim = maximum(radii) + 0.4
    limits!(ax, -lim, lim, -lim, lim)
    return plt
end

@testitem "plot(::HelicalWheel)" begin
    using PDBTools, LazyArtifacts
    import CairoMakie

    seq = "LKKLLKLLKKLLKL"
    n = length(seq)
    fig = CairoMakie.plot(helical_wheel(seq, collect(1.0:n)); colorbar_title = "CN")
    @test fig isa CairoMakie.Figure
    tmpfile = tempname() * ".png"
    CairoMakie.save(tmpfile, fig)
    @test isfile(tmpfile)
    # more than 18 residues (second ring), fixed color limits, no arrow, no color bar
    fig = CairoMakie.plot(helical_wheel(repeat(seq, 2), rand(2n)); normalize = false, clims = (0, 1),
        show_moment = false, colorbar = false, size = (400, 400))
    @test Tuple(fig.scene.viewport[].widths) == (400, 400)
    # all values equal (e.g. zero B-factors), with and without normalization
    @test CairoMakie.plot(helical_wheel(seq, zeros(n))) isa CairoMakie.Figure
    @test CairoMakie.plot(helical_wheel(seq, zeros(n)); normalize = false) isa CairoMakie.Figure
    # labels are white on dark markers and black on light ones
    label_color = Base.get_extension(PDBTools, :MakiePlotting)._label_color
    gradient = CairoMakie.cgrad(:tempo)
    @test label_color(gradient[0.0]) == :black
    @test label_color(gradient[1.0]) == :white
    @test label_color(:black) == :white
    @test label_color(:white) == :black

    # from a selection of atoms
    pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    atoms = select(pdb, at -> resnum(at) <= 18)
    @test CairoMakie.plot(helical_wheel(atoms); colorbar_title = "B-factor") isa CairoMakie.Figure
    # into a position of a figure, and into an axis
    fig = CairoMakie.Figure()
    @test CairoMakie.plot(fig[1, 1], helical_wheel(atoms)) isa CairoMakie.GridLayout
    @test CairoMakie.plot!(CairoMakie.Axis(fig[1, 2]), helical_wheel(atoms)) isa CairoMakie.AbstractPlot
end
