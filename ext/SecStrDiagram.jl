module SecStrDiagram

using TestItems: @testitem
using PDBTools
using PDBTools: Atom, _ss_runs, _ss_transform_contrib, _ss_assign, _SS_HELIX_CODES, _SS_BETA_CODES
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

"""
    find_cmu_fonts()

Look for the Computer Modern Unicode font files (CMU Serif): `cmunrm` (regular), `cmunbx`
(bold) and `cmunti` (italic), as `.ttf` or `.otf`, in the system and user font folders
(Linux, macOS, Windows). Returns a NamedTuple with the paths found (empty if the regular
font is not found).
"""
function find_cmu_fonts()
    dirs = [joinpath(homedir(), ".local", "share", "fonts"),
            joinpath(homedir(), ".fonts"),
            "/usr/local/share/fonts",
            "/usr/share/fonts",
            joinpath(homedir(), "Library", "Fonts"),
            "/Library/Fonts",
            joinpath(get(ENV, "LOCALAPPDATA", ""), "Microsoft", "Windows", "Fonts"),
            joinpath(get(ENV, "WINDIR", "C:\\Windows"), "Fonts")]
    found = Dict{String,String}()
    for d in dirs
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

function PDBTools.ss(pdbfile::AbstractString,
    selection::AbstractString = "protein",
    contrib::Union{Nothing,AbstractVector{<:Real}} = nothing;
    kwargs...,
)
    atoms = if lowercase(splitext(pdbfile)[2]) in (".cif", ".mmcif")
        read_mmcif(pdbfile, selection)
    else
        read_pdb(pdbfile, selection)
    end
    return PDBTools.ss(atoms, contrib; kwargs...)
end

function PDBTools.ss(ats::AbstractVector{<:Atom},
    contrib::Union{Nothing,AbstractVector{<:Real}} = nothing;
    color = nothing,
    normalize::Symbol = :robust,
    range_mode::Symbol = :auto,
    scale::Union{Nothing,Real} = nothing,
    limit::Union{Nothing,Real} = nothing,
    clip_quantile::Real = 0.99,
    signed_tol::Real = 1e-6,
    beta_reduce::Symbol = :mean,
    per_chain::Bool = false,
    chain_new_row::Bool = true,
    colorbar::Bool = true,
    colorbar_label::AbstractString = "",
    output::Union{Nothing,AbstractString} = nothing,
)
    # --- 1) Residues ---
    residues = collect(eachresidue(ats))
    nres = length(residues)
    nres > 0 || throw(ArgumentError("The atom selection has no residues."))

    chain_ids = [string(chain(first(r))) for r in residues]
    resnums = [Int(resnum(first(r))) for r in residues]

    # Start of a segment: first position, change of chain or gap in the numbering
    newseg = [k == 1 ||
              chain_ids[k] != chain_ids[k-1] ||
              resnums[k] - resnums[k-1] != 1 for k in 1:nres]

    # Contiguous blocks of each chain
    chain_blocks = UnitRange{Int}[]
    start = 1
    for k in 2:(nres+1)
        if k == nres + 1 || chain_ids[k] != chain_ids[k-1]
            push!(chain_blocks, start:(k-1))
            start = k
        end
    end
    multichain = length(unique(chain_ids)) > 1

    # --- 2) Per-residue contributions ---
    contrib_raw = if isnothing(contrib)
        betas = [Float64[beta(atom) for atom in r] for r in residues]
        if beta_reduce == :mean
            [sum(b) / length(b) for b in betas]
        elseif beta_reduce == :sum
            # If all atoms of each residue have the same value, beta is probably already a
            # per-residue value; summing would multiply it by the number of atoms.
            if all(b -> all(≈(first(b)), b), betas) && any(b -> length(b) > 1, betas)
                @warn "beta_reduce = :sum, but the beta values are constant within each " *
                      "residue (they look like per-residue values). Consider :mean."
            end
            sum.(betas)
        else
            throw(ArgumentError("beta_reduce must be :mean or :sum, got $beta_reduce"))
        end
    else
        length(contrib) == nres || throw(ArgumentError(
            "Length of contrib ($(length(contrib))) does not match the number of residues ($nres)."
        ))
        Float64.(contrib)
    end

    if per_chain
        for b in chain_blocks
            m = maximum(abs, view(contrib_raw, b))
            iszero(m) || (contrib_raw[b] ./= m)
        end
    end

    # --- 3) Secondary structure (DSSP) ---
    ssvec = _ss_assign(ats, residues)

    # Helices shorter than 4 residues and strands shorter than 3 become loops
    for (i, j) in _ss_runs([c in _SS_HELIX_CODES for c in ssvec], newseg)
        (j - i + 1) < 4 && (ssvec[i:j] .= "C")
    end
    for (i, j) in _ss_runs([c in _SS_BETA_CODES for c in ssvec], newseg)
        (j - i + 1) < 3 && (ssvec[i:j] .= "C")
    end

    is_helix = [c in _SS_HELIX_CODES for c in ssvec]
    is_beta = [c in _SS_BETA_CODES for c in ssvec]

    # --- 4) Colors ---
    contrib_plot, crange, s = _ss_transform_contrib(contrib_raw;
        mode = normalize, range_mode, scale, limit, clip_quantile, signed_tol)
    signed = crange[1] < 0
    cmap = cgrad(isnothing(color) ? (signed ? :balance : :tempo) : color)
    # Values out of the range are painted with the color of the extreme of the color map.

    # --- 5) Geometry (global coordinates: residue k occupies x ∈ [k, k+1], y ∈ [0, 3]) ---
    # Layers: 1 = loops, 2 = helix (odd position in the segment, behind),
    #         3 = helix (even position in the segment, in front), 4 = β strands
    y_bottom, y_top = 0.0, 3.0
    shapes = Tuple{Int,Int,Vector{Point2f}}[] # (layer, residue, vertices)

    for k in 1:nres
        if !is_helix[k] && !is_beta[k]
            push!(shapes, (1, k, Point2f[(k, y_bottom + 1), (k + 1, y_bottom + 1),
                                         (k + 1, y_top - 1), (k, y_top - 1)]))
        end
    end

    for (i, j) in _ss_runs(is_helix, newseg)
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

    for (i, j) in _ss_runs(is_beta, newseg)
        for k in i:(j-1)
            push!(shapes, (4, k, Point2f[(k, y_bottom + 0.8), (k + 1, y_bottom + 0.8),
                                         (k + 1, y_top - 0.8), (k, y_top - 0.8)]))
        end
        y_mid = (y_bottom + y_top) / 2
        push!(shapes, (4, j, Point2f[(j, y_bottom), (j + 1, y_mid), (j, y_top)]))
    end

    sort!(shapes; by = first)

    # --- 6) Rows of the figure ---
    blocks = chain_new_row ? chain_blocks : [1:nres]
    rows = UnitRange{Int}[]
    for b in blocks, lo in first(b):RESIDUES_PER_ROW:last(b)
        push!(rows, lo:min(lo + RESIDUES_PER_ROW - 1, last(b)))
    end
    ncols = maximum(length, rows)

    labels = PDBTools.residue_ticks(ats; oneletter = true)[2]
    length(labels) == nres || error("Residue labels ($(length(labels))) do not match residues ($nres).")

    # --- 7) Font (Computer Modern, if found) ---
    fonts = find_cmu_fonts()
    isempty(fonts) && @warn "Computer Modern (cmunrm.ttf) not found; using the default " *
        "Makie font. Install the CMU Serif fonts to use Computer Modern." maxlog = 1
    fig_kwargs = isempty(fonts) ? (;) : (; fonts = fonts)

    # The size of the figure is adjusted to the content at the end (resize_to_layout!),
    # from the fixed width per residue and the fixed height per row.
    f = Figure(; fontsize = FONTSIZE, figure_padding = FIGURE_PADDING,
        backgroundcolor = :transparent, fig_kwargs...)

    # 0.5 residue of space on each side for the helices at the edges
    xspan = ncols + 1.0

    # --- 8) Plot, row by row ---
    for (r, idx) in enumerate(rows)
        off = first(idx) - 1
        n = length(idx)
        ax = Axis(f[r, 1];
            width = xspan * PT_PER_RESIDUE,
            height = ROW_HEIGHT,
            xticks = (collect(1:n) .+ 0.5, labels[idx]),
            xticklabelrotation = pi / 2,
            ylabel = (multichain && chain_new_row) ? "Chain $(chain_ids[first(idx)])" : "",
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
                color = contrib_plot[k],
                colormap = cmap,
                colorrange = crange,
                strokecolor = :black,
                strokewidth = STROKEWIDTH)
        end
    end
    rowgap!(f.layout, ROW_GAP)

    # --- 9) Horizontal color bar below the rows (ticks in original units) ---
    if colorbar
        tick_pos = collect(range(crange[1], crange[2]; length = 5))
        tick_val = normalize == :asinh ? sinh.(tick_pos) .* s : tick_pos
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
            label = colorbar_label,
            ticks = (tick_pos, tick_lab))
        rowgap!(f.layout, length(rows), 3 * ROW_GAP)
    end

    resize_to_layout!(f)

    # Always SVG (vector); 1 unit = 1 pt
    if !isnothing(output)
        svgfile = splitext(output)[1] * ".svg"
        svgfile != output && @info "Saving as SVG: $svgfile"
        CairoMakie.save(svgfile, f; pt_per_unit = 1)
    end
    return f
end

@testitem "ss diagram" begin
    using PDBTools, LazyArtifacts
    import CairoMakie
    atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    nres = length(eachresidue(atoms))

    f = ss(atoms)
    @test f isa CairoMakie.Figure

    tmpfile = tempname() * ".png"
    f = ss(atoms, collect(range(-1.0, 1.0; length = nres)); normalize = :asinh, output = tmpfile)
    @test isfile(splitext(tmpfile)[1] * ".svg")

    f = ss(PDBTools.test_file("structure.pdb"), "protein"; beta_reduce = :sum, per_chain = true)
    @test f isa CairoMakie.Figure

    @test_throws ArgumentError ss(atoms, ones(nres - 1))
    @test_throws ArgumentError ss(atoms; beta_reduce = :median)
end

end # module SecStrDiagram
