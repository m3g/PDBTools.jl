#
# Secondary structure diagrams: plot(::SecondaryStructureDiagram)
#
# --- Fixed layout of the figure (units = typographic points; 72 pt = 1 inch) ---
# Rows of 50 residues × 13.5 pt = 675 pt ≈ 24 cm, with a 12 pt font. Saved with
# `pt_per_unit = 2/3`, the width is that of a full page (≈ 16 cm), with an 8 pt font.
const RESIDUES_PER_ROW = 50    # residues per row
const PT_PER_RESIDUE = 13.5    # width of each residue (50 × 13.5 = 675 pt)
const ROW_HEIGHT = 60.0        # height of the drawing area of each row (without labels)
const ROW_GAP = 6.0            # space between rows
const FONTSIZE = 12            # font size, in pt
const STROKEWIDTH = 0.75       # outline of the elements, in pt
const COLORBAR_WIDTH = 300.0   # length of the horizontal color bar, in pt
const COLORBAR_HEIGHT = 9.0    # thickness of the color bar, in pt
const FIGURE_PADDING = 6.0     # outer margin, in pt
const NO_VALUE_COLOR = :grey85 # color of the residues when there are no values

const HELIX_CODES = ("G", "H", "I", "P")
const STRAND_CODES = ("E", "B")

# Contiguous segments (i, j) where mask is true. A segment is broken when
# segment_start[k] is true (start of a new chain or a gap in the numbering), so that
# a helix or strand never crosses the boundary between two chains.
function _structure_runs(mask::AbstractVector{Bool}, segment_start::AbstractVector{Bool})
    length(mask) == length(segment_start) ||
        throw(ArgumentError("mask and segment_start must have the same length"))
    runs = Tuple{Int,Int}[]
    n = length(mask)
    i = 1
    while i <= n
        if mask[i]
            j = i
            while j < n && mask[j+1] && !segment_start[j+1]
                j += 1
            end
            push!(runs, (i, j))
            i = j + 1
        else
            i += 1
        end
    end
    return runs
end

# Start of a segment: first residue, change of chain or gap in the numbering
function _segment_starts(chains::AbstractVector, residue_numbers::AbstractVector{<:Integer})
    return [k == 1 || chains[k] != chains[k-1] || residue_numbers[k] - residue_numbers[k-1] != 1
            for k in eachindex(chains)]
end

# Sample quantile with linear interpolation (same definition as the default of Statistics.quantile)
function _sample_quantile(v::AbstractVector{<:Real}, p::Real)
    0 <= p <= 1 || throw(ArgumentError("quantile must be in [0, 1], got $p"))
    s = sort(v)
    n = length(s)
    n == 1 && return float(s[1])
    h = (n - 1) * p + 1
    lo = floor(Int, h)
    hi = min(lo + 1, n)
    return s[lo] + (h - lo) * (s[hi] - s[lo])
end

# Transform the values for the color mapping. Returns (y, crange, s): transformed values,
# range of the color scale (in transformed units) and scale s (relevant only for :asinh).
function _color_values(x::AbstractVector{<:Real};
    color_scale::Symbol = :robust,
    range_mode::Symbol = :auto,
    asinh_scale::Union{Nothing,Real} = nothing,
    color_limit::Union{Nothing,Real} = nothing,
    clip_quantile::Real = 0.99,
    signed_tolerance::Real = 1e-6,
)
    xf = float.(x)
    ref = maximum(abs, xf)

    # Limit in original units
    lim0 = if !isnothing(color_limit)
        float(color_limit)
    elseif color_scale == :robust
        q = _sample_quantile(abs.(xf), clip_quantile)
        iszero(q) ? ref : q # many zeros: use the maximum
    elseif color_scale in (:linear, :asinh)
        ref
    else
        throw(ArgumentError("color_scale must be :robust, :linear or :asinh, got $color_scale"))
    end
    lim0 = (iszero(lim0) || !isfinite(lim0)) ? 1.0 : lim0 # all values are zero
    lim0 > 0 || throw(ArgumentError("color_limit must be positive, got $lim0"))

    # Type of range
    signed = if range_mode == :auto
        minimum(xf) < -signed_tolerance * ref
    elseif range_mode == :symmetric
        true
    elseif range_mode == :positive
        false
    else
        throw(ArgumentError("range_mode must be :auto, :symmetric or :positive, got $range_mode"))
    end

    # Transformation
    s = 1.0
    if color_scale == :asinh
        s = isnothing(asinh_scale) ? lim0 / 20 : float(asinh_scale)
        s > 0 || throw(ArgumentError("asinh_scale must be positive, got $s"))
        y = asinh.(xf ./ s)
        lim = asinh(lim0 / s)
    else
        y = xf
        lim = lim0
    end

    crange = signed ? (-lim, lim) : (zero(lim), lim)
    return y, crange, s
end

"""
    plot(diagram::SecondaryStructureDiagram; kwargs...)

Plot a secondary structure diagram, as computed by [`secondary_structure_diagram`](@ref): the secondary
structure of each residue along the sequence (helices as ribbons, β strands as arrows and
loops as thin bars), with each residue colored by its value.

!!! note
    This function requires a Makie backend, e.g. `using CairoMakie`.

Returns a `Makie.Figure`. Helices shorter than 4 residues and strands shorter than 3
residues are drawn as loops.

# Layout

The figure has a fixed layout, so that all diagrams come out identical: 50 residues per
row, 13.5 pt per residue (≈ 24 cm wide rows) and a 12 pt font. The height grows with the
number of rows. To keep the sizes in typographic points, save it as SVG with CairoMakie:
`CairoMakie.save("diagram.svg", fig; pt_per_unit = 1)`. With `pt_per_unit = 2/3`, the
figure has the width of a full page (≈ 16 cm), with an 8 pt font.

# Keyword arguments
- `colormap = nothing`: color map. If `nothing`, uses `:balance` (diverging, zero = white)
  when the color range is symmetric and `:tempo` (sequential) when it is `(0, limit)`.
- `colorbar = true`, `colorbar_title = ""`: horizontal color bar below the rows (not drawn
  when there are no values).
- `color_scale = :robust`: how values are mapped to colors:
  - `:robust`: linear, with the limit at the quantile `clip_quantile` of the absolute values;
  - `:linear`: linear, with the limit at the maximum absolute value;
  - `:asinh`: `asinh(value / asinh_scale)`, linear for small values and logarithmic for
    large ones, for values that span several orders of magnitude.
- `range_mode = :auto`: `:symmetric` (`(-limit, limit)`, zero at the center of the color
  map), `:positive` (`(0, limit)`) or `:auto` (symmetric if there are negative values).
- `color_limit = nothing`, `asinh_scale = nothing`: limit of the color scale and scale of
  `:asinh`, in the units of the values. Use the same values to compare several diagrams.
- `clip_quantile = 0.99`: quantile used by `:robust`.
- `signed_tolerance = 1e-6`: values smaller than `-signed_tolerance * maximum(abs, values)`
  are considered negative by `range_mode = :auto`.
- `normalize_per_chain = false`: divide the values of each chain by the maximum absolute
  value of that chain.
- `chain_new_row = true`: each chain starts in a new row, labeled on the left.
- `background_color = :white`: background of the figure (e.g. `:transparent`).
- `font = nothing`: font of the labels (default: New Computer Modern, distributed with Makie).
- `size = nothing`: size of the figure, as `(width, height)`, in pt. By default, the
  size is fixed by the content. With a given size, the rows are stretched to fill the
  figure (the font size is not changed).

Values beyond the color limits take the color of the extreme of the color map.

# Example

```julia-repl
julia> using PDBTools, LazyArtifacts, CairoMakie

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> fig = plot(secondary_structure_diagram(atoms); colorbar_title = "B-factor")

julia> CairoMakie.save("diagram.svg", fig; pt_per_unit = 1)
```

"""
function Makie.plot(d::SecondaryStructureDiagram;
    size = nothing,
    font = nothing,
    background_color = :white,
    kwargs...,
)
    # By default, the size of the figure is adjusted to the content at the end
    # (resize_to_layout!), from the fixed width per residue and the fixed height per row.
    # With a given `size`, the rows are stretched to fill the figure.
    f = Figure(; fontsize = FONTSIZE, figure_padding = FIGURE_PADDING,
        backgroundcolor = background_color, fonts = _fonts(font),
        (isnothing(size) ? (;) : (; size = _check_size(size)))...)
    _draw_ss!(f.layout, d; fixed = isnothing(size), background_color, kwargs...)
    isnothing(size) && resize_to_layout!(f)
    return f
end

function Makie.plot(gp::GridPos, d::SecondaryStructureDiagram; background_color = :white, kwargs...)
    gl = GridLayout(gp)
    _draw_ss!(gl, d; fixed = false, background_color, kwargs...)
    return gl
end

function Makie.plot!(::Axis, ::SecondaryStructureDiagram; kwargs...)
    throw(ArgumentError("""\n
        A secondary structure diagram has one axis per row, and cannot be drawn into an existing Axis.
        Use plot(fig[i, j], diagram) to draw it into a position of a figure.

    """))
end

# Draw the diagram into the layout `gl`: rows of residues, and the color bar below. With
# `fixed = true`, the axes have the fixed size of the layout; otherwise they fill `gl`.
function _draw_ss!(gl::GridLayout, d::SecondaryStructureDiagram;
    fixed::Bool,
    colormap = nothing,
    colorbar::Bool = true,
    colorbar_title::AbstractString = "",
    color_scale::Symbol = :robust,
    range_mode::Symbol = :auto,
    color_limit::Union{Nothing,Real} = nothing,
    asinh_scale::Union{Nothing,Real} = nothing,
    clip_quantile::Real = 0.99,
    signed_tolerance::Real = 1e-6,
    normalize_per_chain::Bool = false,
    chain_new_row::Bool = true,
    background_color = :white,
)
    ss = d.ss
    nres = length(ss)
    nres > 0 || throw(ArgumentError("The secondary structure data is empty."))

    chains = [string(s.chain) for s in ss]
    residue_numbers = [Int(s.resnum) for s in ss]
    labels = [string(PDBTools.oneletter(s.resname), s.resnum) for s in ss]
    segment_start = _segment_starts(chains, residue_numbers)
    multichain = length(unique(chains)) > 1

    # Secondary structure; short helices and strands are drawn as loops
    codes = [string(s.sscode) for s in ss]
    for (i, j) in _structure_runs([c in HELIX_CODES for c in codes], segment_start)
        (j - i + 1) < 4 && (codes[i:j] .= "C")
    end
    for (i, j) in _structure_runs([c in STRAND_CODES for c in codes], segment_start)
        (j - i + 1) < 3 && (codes[i:j] .= "C")
    end
    is_helix = [c in HELIX_CODES for c in codes]
    is_strand = [c in STRAND_CODES for c in codes]

    # Contiguous blocks of each chain
    chain_blocks = UnitRange{Int}[]
    start = 1
    for k in 2:(nres+1)
        if k == nres + 1 || chains[k] != chains[k-1]
            push!(chain_blocks, start:(k-1))
            start = k
        end
    end

    # --- Colors ---
    has_values = !isnothing(d.values)
    if has_values
        values = copy(d.values)
        if normalize_per_chain
            for b in chain_blocks
                m = maximum(abs, view(values, b))
                iszero(m) || (values[b] ./= m)
            end
        end
        color_values, crange, s = _color_values(values;
            color_scale, range_mode, asinh_scale, color_limit, clip_quantile, signed_tolerance)
        cmap = cgrad(isnothing(colormap) ? (crange[1] < 0 ? :balance : :tempo) : colormap)
    end
    residue_color(k) = has_values ?
        (; color = color_values[k], colormap = cmap, colorrange = crange) :
        (; color = NO_VALUE_COLOR)

    # --- Geometry (global coordinates: residue k occupies x ∈ [k, k+1], y ∈ [0, 3]) ---
    # Layers: 1 = loops, 2 = helix (odd position in the segment, behind),
    #         3 = helix (even position in the segment, in front), 4 = β strands
    y_bottom, y_top = 0.0, 3.0
    shapes = Tuple{Int,Int,Vector{Point2f}}[] # (layer, residue, vertices)

    for k in 1:nres
        if !is_helix[k] && !is_strand[k]
            push!(shapes, (1, k, Point2f[(k, y_bottom + 1), (k + 1, y_bottom + 1),
                                         (k + 1, y_top - 1), (k, y_top - 1)]))
        end
    end

    for (i, j) in _structure_runs(is_helix, segment_start)
        for k in i:j
            x1, x2, x3 = k - 0.5, k + 0.5, k + 1.5
            idx_in_run = k - i
            is_first = (k == i)
            is_last = (k == j)
            yb, yt = if is_first
                (y_bottom + 1.0, y_top - 0.5)
            elseif is_last && iseven(idx_in_run)
                (y_bottom + 0.5, y_top - 1.0)
            elseif is_last
                (y_bottom + 1.0, y_top - 0.5)
            else
                (y_bottom + 0.5, y_top - 0.5)
            end
            pts = if iseven(idx_in_run)
                Point2f[(x1, yb), (x2, yb), (x3, yt), (x2, yt)] # tilted to the right
            else
                Point2f[(x1, yt), (x2, yt), (x3, yb), (x2, yb)] # tilted to the left
            end
            push!(shapes, (iseven(idx_in_run) ? 3 : 2, k, pts))
        end
    end

    for (i, j) in _structure_runs(is_strand, segment_start)
        for k in i:(j-1)
            push!(shapes, (4, k, Point2f[(k, y_bottom + 0.8), (k + 1, y_bottom + 0.8),
                                         (k + 1, y_top - 0.8), (k, y_top - 0.8)]))
        end
        y_mid = (y_bottom + y_top) / 2
        push!(shapes, (4, j, Point2f[(j, y_bottom), (j + 1, y_mid), (j, y_top)]))
    end

    sort!(shapes; by = first)

    # --- Rows of the figure ---
    blocks = chain_new_row ? chain_blocks : [1:nres]
    rows = UnitRange{Int}[]
    for b in blocks, lo in first(b):RESIDUES_PER_ROW:last(b)
        push!(rows, lo:min(lo + RESIDUES_PER_ROW - 1, last(b)))
    end
    ncols = maximum(length, rows)

    # 0.5 residue of space on each side for the helices at the edges
    xspan = ncols + 1.0

    # --- Plot, row by row ---
    for (r, idx) in enumerate(rows)
        off = first(idx) - 1
        n = length(idx)
        ax = Axis(gl[r, 1];
            (fixed ? (; width = xspan * PT_PER_RESIDUE, height = ROW_HEIGHT) : (;))...,
            xticks = (collect(1:n) .+ 0.5, labels[idx]),
            xticklabelrotation = pi / 2,
            ylabel = (multichain && chain_new_row) ? "Chain $(chains[first(idx)])" : "",
            backgroundcolor = background_color)
        ax.xgridvisible = false
        ax.ygridvisible = false
        ax.yticklabelsvisible = false
        ax.yticksvisible = false
        ax.xticksvisible = false
        hidespines!(ax)
        # same scale in all rows (short rows are aligned to the left)
        limits!(ax, 0.5, 0.5 + xspan, y_bottom, y_top)

        for (_, k, pts) in shapes
            first(idx) <= k <= last(idx) || continue
            poly!(ax, [p - Point2f(off, 0) for p in pts];
                residue_color(k)...,
                strokecolor = :black,
                strokewidth = STROKEWIDTH)
        end
    end
    rowgap!(gl, ROW_GAP)

    # --- Horizontal color bar below the rows (ticks in the units of the values) ---
    if colorbar && has_values
        tick_pos = collect(range(crange[1], crange[2]; length = 5))
        tick_val = color_scale == :asinh ? sinh.(tick_pos) .* s : tick_pos
        # zero the floating point residue in the central tick (e.g. 1e-17)
        tol = 1e-9 * maximum(abs, tick_val)
        tick_lab = [abs(v) <= tol ? "0" : string(round(v; sigdigits = 2)) for v in tick_val]
        Colorbar(gl[length(rows)+1, 1];
            vertical = false,
            flipaxis = false,
            width = COLORBAR_WIDTH,
            height = COLORBAR_HEIGHT,
            tellwidth = false,
            colormap = cmap,
            limits = crange,
            label = colorbar_title,
            ticks = (tick_pos, tick_lab))
        rowgap!(gl, length(rows), 3 * ROW_GAP)
    end
    return gl
end

@testitem "plot(::SecondaryStructureDiagram) helpers" begin
    using PDBTools
    import CairoMakie
    ext = Base.get_extension(PDBTools, :MakiePlotting)

    runs = ext._structure_runs
    @test runs(Bool[0, 1, 1, 0, 1], falses(5)) == [(2, 3), (5, 5)]
    # a new segment breaks a run
    @test runs(Bool[1, 1, 1, 1], Bool[1, 0, 1, 0]) == [(1, 2), (3, 4)]
    @test isempty(runs(falses(3), falses(3)))
    @test_throws ArgumentError runs(Bool[1, 0], Bool[1])

    @test ext._segment_starts(["A", "A", "A", "B"], [1, 2, 4, 5]) == [true, false, true, true]

    @test ext._sample_quantile([1, 2, 3, 4, 5], 0.5) == 3.0
    @test ext._sample_quantile([1, 2, 3, 4], 0.5) == 2.5
    @test ext._sample_quantile([7], 0.9) == 7.0
    @test_throws ArgumentError ext._sample_quantile([1, 2], 1.5)

    cv = ext._color_values
    x = collect(0.0:10.0)
    y, crange, s = cv(x; color_scale = :linear)
    @test y == x
    @test crange == (0.0, 10.0)
    # negative values give a symmetric range; numerical noise does not
    @test cv([-5.0, 10.0]; color_scale = :linear)[2] == (-10.0, 10.0)
    @test cv([-1e-12, 10.0]; color_scale = :linear)[2] == (0.0, 10.0)
    @test cv(x; color_scale = :linear, range_mode = :symmetric)[2] == (-10.0, 10.0)
    @test cv([-5.0, 10.0]; color_scale = :linear, range_mode = :positive)[2] == (0.0, 10.0)
    # robust limit is the quantile of |x|
    @test cv(x; color_scale = :robust, clip_quantile = 0.5)[2] == (0.0, 5.0)
    # fixed limit; all-zero data do not break the scale
    @test cv(x; color_limit = 2.0)[2] == (0.0, 2.0)
    @test cv(zeros(4))[2] == (0.0, 1.0)
    # asinh
    y, crange, s = cv(x; color_scale = :asinh, asinh_scale = 1.0)
    @test y ≈ asinh.(x)
    @test crange[2] ≈ asinh(10.0)
    @test s == 1.0
    @test cv(x; color_scale = :asinh)[3] ≈ 10.0 / 20
    @test_throws ArgumentError cv(x; color_scale = :log)
    @test_throws ArgumentError cv(x; range_mode = :negative)
    @test_throws ArgumentError cv(x; color_limit = -1.0)
    @test_throws ArgumentError cv(x; color_scale = :asinh, asinh_scale = -1.0)
end

@testitem "plot(::SecondaryStructureDiagram)" begin
    using PDBTools, LazyArtifacts
    import CairoMakie
    atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    nres = length(eachresidue(atoms))

    # default property (B-factors)
    d = secondary_structure_diagram(atoms)
    fig = CairoMakie.plot(d; colorbar_title = "B-factor")
    @test fig isa CairoMakie.Figure
    svgfile = tempname() * ".svg"
    CairoMakie.save(svgfile, fig; pt_per_unit = 1)
    @test isfile(svgfile)

    # signed values, color options, layouts
    ds = secondary_structure_diagram(atoms, collect(range(-1.0, 1.0; length = nres)))
    @test CairoMakie.plot(ds; color_scale = :asinh) isa CairoMakie.Figure
    @test CairoMakie.plot(ds; color_scale = :linear, color_limit = 0.5,
        colormap = :viridis, background_color = :transparent, font = "Arial") isa CairoMakie.Figure
    @test CairoMakie.plot(ds; normalize_per_chain = true) isa CairoMakie.Figure
    @test CairoMakie.plot(ds; chain_new_row = false, colorbar = false) isa CairoMakie.Figure
    # fixed size by default; a given size is respected
    fig = CairoMakie.plot(ds; size = (800, 600))
    @test Tuple(fig.scene.viewport[].widths) == (800, 600)
    @test Tuple(CairoMakie.plot(ds).scene.viewport[].widths) != (800, 600)
    @test_throws ArgumentError CairoMakie.plot(ds; size = (0, 600))

    # only the secondary structure
    @test CairoMakie.plot(secondary_structure_diagram(dssp_run(atoms))) isa CairoMakie.Figure

    # into a position of a figure; not into an Axis
    fig = CairoMakie.Figure(; size = (900, 600))
    @test CairoMakie.plot(fig[1, 1], d) isa CairoMakie.GridLayout
    @test_throws ArgumentError CairoMakie.plot!(CairoMakie.Axis(fig[2, 1]), d)

    @test_throws ArgumentError CairoMakie.plot(ds; color_scale = :log)
    @test_throws ArgumentError CairoMakie.plot(ds; range_mode = :negative)
end
