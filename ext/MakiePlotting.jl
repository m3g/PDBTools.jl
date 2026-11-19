module MakiePlotting

using TestItems: @testitem
using PDBTools
using PDBTools: SecondaryStructureDiagram, _structure_runs, _segment_starts, _color_values,
    _HELIX_CODES, _STRAND_CODES
import CairoMakie: plot
using CairoMakie: CairoMakie, Figure, Axis, Colorbar, Point2f, poly!, limits!, hidespines!,
    rowgap!, resize_to_layout!, cgrad

# --- Fixed layout of the figure (units = typographic points; 72 pt = 1 inch) ---
# Usable width of a full page ≈ 17 cm ≈ 482 pt.
const RESIDUES_PER_ROW = 50    # residues per row
const PT_PER_RESIDUE = 9.0     # width of each residue (50 × 9 = 450 pt)
const ROW_HEIGHT = 40.0        # height of the drawing area of each row (without labels)
const ROW_GAP = 4.0            # space between rows
const FONTSIZE = 8             # font size, in pt
const STROKEWIDTH = 0.5        # outline of the elements, in pt
const COLORBAR_WIDTH = 200.0   # length of the horizontal color bar, in pt
const COLORBAR_HEIGHT = 6.0    # thickness of the color bar, in pt
const FIGURE_PADDING = 4.0     # outer margin, in pt

# System and user font folders (Linux, macOS, Windows)
_font_directories() = [joinpath(homedir(), ".local", "share", "fonts"),
                       joinpath(homedir(), ".fonts"),
                       "/usr/local/share/fonts",
                       "/usr/share/fonts",
                       joinpath(homedir(), "Library", "Fonts"),
                       "/Library/Fonts",
                       joinpath(get(ENV, "LOCALAPPDATA", ""), "Microsoft", "Windows", "Fonts"),
                       joinpath(get(ENV, "WINDIR", "C:\\Windows"), "Fonts")]

"""
    _find_computer_modern_fonts(directories = _font_directories())

Look for the Computer Modern Unicode font files (CMU Serif): `cmunrm` (regular), `cmunbx`
(bold) and `cmunti` (italic), as `.ttf` or `.otf`, in `directories` (by default, the
system and user font folders). Returns a NamedTuple with the paths found (empty if the
regular font is not found).
"""
function _find_computer_modern_fonts(directories = _font_directories())
    found = Dict{String,String}()
    for d in directories
        isdir(d) || continue
        for (root, _, files) in walkdir(d), file in files
            name, ext = splitext(lowercase(file))
            ext in (".ttf", ".otf") || continue
            name in ("cmunrm", "cmunbx", "cmunti") || continue
            haskey(found, name) || (found[name] = joinpath(root, file))
        end
        haskey(found, "cmunrm") && break
    end
    haskey(found, "cmunrm") || return (;)
    fonts = (; regular = found["cmunrm"])
    haskey(found, "cmunbx") && (fonts = merge(fonts, (; bold = found["cmunbx"])))
    haskey(found, "cmunti") && (fonts = merge(fonts, (; italic = found["cmunti"])))
    return fonts
end

# Plotting of a SecondaryStructureDiagram. The keyword arguments are documented in the
# docstring of PDBTools.SecondaryStructureDiagram.
function plot(diagram::SecondaryStructureDiagram;
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
)
    chains = diagram.chains
    nres = length(chains)
    segment_start = _segment_starts(chains, diagram.residue_numbers)
    multichain = length(unique(chains)) > 1

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
    values = copy(diagram.values)
    if normalize_per_chain
        for b in chain_blocks
            m = maximum(abs, view(values, b))
            iszero(m) || (values[b] ./= m)
        end
    end
    color_values, crange, s = _color_values(values;
        color_scale, range_mode, asinh_scale, color_limit, clip_quantile, signed_tolerance)
    signed = crange[1] < 0
    cmap = cgrad(isnothing(colormap) ? (signed ? :balance : :tempo) : colormap)

    is_helix = [c in _HELIX_CODES for c in diagram.secondary_structure]
    is_strand = [c in _STRAND_CODES for c in diagram.secondary_structure]

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
    labels = diagram.residue_labels

    # --- Font (Computer Modern, if found) ---
    fonts = _find_computer_modern_fonts()
    isempty(fonts) && @warn "Computer Modern (cmunrm.ttf) not found; using the default " *
        "Makie font. Install the CMU Serif fonts to use Computer Modern." maxlog = 1
    fig_kwargs = isempty(fonts) ? (;) : (; fonts = fonts)

    # The size of the figure is adjusted to the content at the end (resize_to_layout!),
    # from the fixed width per residue and the fixed height per row.
    f = Figure(; fontsize = FONTSIZE, figure_padding = FIGURE_PADDING,
        backgroundcolor = :transparent, fig_kwargs...)

    # 0.5 residue of space on each side for the helices at the edges
    xspan = ncols + 1.0

    # --- Plot, row by row ---
    for (r, idx) in enumerate(rows)
        off = first(idx) - 1
        n = length(idx)
        ax = Axis(f[r, 1];
            width = xspan * PT_PER_RESIDUE,
            height = ROW_HEIGHT,
            xticks = (collect(1:n) .+ 0.5, labels[idx]),
            xticklabelrotation = pi / 2,
            ylabel = (multichain && chain_new_row) ? "Chain $(chains[first(idx)])" : "",
            backgroundcolor = :transparent)
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
                color = color_values[k],
                colormap = cmap,
                colorrange = crange,
                strokecolor = :black,
                strokewidth = STROKEWIDTH)
        end
    end
    rowgap!(f.layout, ROW_GAP)

    # --- Horizontal color bar below the rows (ticks in the units of the values) ---
    if colorbar
        tick_pos = collect(range(crange[1], crange[2]; length = 5))
        tick_val = color_scale == :asinh ? sinh.(tick_pos) .* s : tick_pos
        # zero the floating point residue in the central tick (e.g. 1e-17)
        tol = 1e-9 * maximum(abs, tick_val)
        tick_lab = [abs(v) <= tol ? "0" : string(round(v; sigdigits = 2)) for v in tick_val]
        Colorbar(f[length(rows)+1, 1];
            vertical = false,
            flipaxis = false,
            width = COLORBAR_WIDTH,
            height = COLORBAR_HEIGHT,
            tellwidth = false,
            colormap = cmap,
            limits = crange,
            label = colorbar_title,
            ticks = (tick_pos, tick_lab))
        rowgap!(f.layout, length(rows), 3 * ROW_GAP)
    end

    resize_to_layout!(f)
    return f
end

@testitem "plot(::SecondaryStructureDiagram)" begin
    using PDBTools, LazyArtifacts
    import CairoMakie
    atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    nres = length(eachresidue(atoms))

    fig = CairoMakie.plot(secondary_structure_diagram(atoms); colorbar_title = "B-factor")
    @test fig isa CairoMakie.Figure
    svgfile = tempname() * ".svg"
    CairoMakie.save(svgfile, fig; pt_per_unit = 1)
    @test isfile(svgfile)

    # signed values, asinh scale, per-chain normalization, one row, no colorbar
    diagram = secondary_structure_diagram(atoms, collect(range(-1.0, 1.0; length = nres)))
    @test CairoMakie.plot(diagram; color_scale = :asinh) isa CairoMakie.Figure
    @test CairoMakie.plot(diagram; color_scale = :linear, color_limit = 0.5, colormap = :viridis) isa CairoMakie.Figure
    @test CairoMakie.plot(diagram; normalize_per_chain = true) isa CairoMakie.Figure
    @test CairoMakie.plot(diagram; chain_new_row = false, colorbar = false) isa CairoMakie.Figure

    @test_throws ArgumentError CairoMakie.plot(diagram; color_scale = :log)
    @test_throws ArgumentError CairoMakie.plot(diagram; range_mode = :negative)
end

@testitem "Computer Modern fonts" begin
    using PDBTools
    import CairoMakie
    ext = Base.get_extension(PDBTools, :MakiePlotting)
    dir = mktempdir()
    @test isempty(ext._find_computer_modern_fonts([dir]))
    for name in ("cmunrm.ttf", "cmunbx.ttf", "cmunti.otf")
        touch(joinpath(dir, name))
    end
    fonts = ext._find_computer_modern_fonts([dir, joinpath(dir, "does_not_exist")])
    @test fonts.regular == joinpath(dir, "cmunrm.ttf")
    @test fonts.bold == joinpath(dir, "cmunbx.ttf")
    @test fonts.italic == joinpath(dir, "cmunti.otf")
end

end # module MakiePlotting
