export secondary_structure_diagram, SecondaryStructureDiagram

const _HELIX_CODES = ("G", "H", "I", "P")
const _STRAND_CODES = ("E", "B")

# Mean of f over the elements of x (Statistics is not a dependency)
_residue_mean(f, x) = sum(f, x) / length(x)

"""
    SecondaryStructureDiagram

Result of [`secondary_structure_diagram`](@ref): the secondary structure and one value
per residue.

# Fields
- `residue_labels::Vector{String}`: one-letter code and number of each residue.
- `chains::Vector{String}`: chain of each residue.
- `residue_numbers::Vector{Int}`: residue numbers.
- `secondary_structure::Vector{String}`: DSSP code of each residue, after converting
  helices shorter than 4 residues and strands shorter than 3 residues to loops (`"C"`).
- `values::Vector{Float64}`: value of each residue.

# Plotting

    plot(diagram::SecondaryStructureDiagram; kwargs...)

Draws the diagram: helices as ribbons, β strands as arrows and loops as thin bars, with
each residue colored by its value. Requires `using CairoMakie`, and returns a
`CairoMakie.Figure`, which can be saved with `save`. If `Plots` is also loaded, call
`CairoMakie.plot(diagram)`.

The figure has a fixed layout, so that all diagrams come out identical and ready for
publication: full page width (≈ 17 cm), 50 residues per row, 9 pt per residue and an
8 pt Computer Modern font (if the CMU Serif font files are installed). The height grows
with the number of rows. To keep the sizes in typographic points, save as SVG with
`save("diagram.svg", fig; pt_per_unit = 1)`.

Keyword arguments:
- `colormap = nothing`: color map. If `nothing`, uses `:balance` (diverging, zero = white)
  when the color range is symmetric and `:tempo` (sequential) when it is `(0, limit)`.
- `colorbar = true`, `colorbar_title = ""`: horizontal color bar below the rows.
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

Values beyond the color limits take the color of the extreme of the color map.

# Example

```julia-repl
julia> using PDBTools, LazyArtifacts, CairoMakie

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> diagram = secondary_structure_diagram(atoms);

julia> fig = plot(diagram; colorbar_title = "B-factor");

julia> save("diagram.svg", fig; pt_per_unit = 1)
```
"""
struct SecondaryStructureDiagram
    residue_labels::Vector{String}
    chains::Vector{String}
    residue_numbers::Vector{Int}
    secondary_structure::Vector{String}
    values::Vector{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", diagram::SecondaryStructureDiagram)
    n = length(diagram.values)
    nchains = length(unique(diagram.chains))
    nhelix = count(c -> c in _HELIX_CODES, diagram.secondary_structure)
    nstrand = count(c -> c in _STRAND_CODES, diagram.secondary_structure)
    vmin, vmax = extrema(diagram.values)
    print(io, chomp("""
    SecondaryStructureDiagram with $n residues in $nchains chain(s)
        helix: $nhelix, strand: $nstrand, loop: $(n - nhelix - nstrand)
        values from $(round(vmin; sigdigits = 4)) to $(round(vmax; sigdigits = 4))
    """))
end

"""
    _structure_runs(mask, segment_start)

Return the contiguous segments `(i, j)` where `mask` is true. A segment is broken
when `segment_start[k]` is true (start of a new chain or a gap in the numbering), so
that a helix or strand never crosses the boundary between two chains.
"""
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

"""
    _color_values(x; color_scale = :robust, range_mode = :auto, asinh_scale = nothing,
                  color_limit = nothing, clip_quantile = 0.99, signed_tolerance = 1e-6)

Transform the per-residue values for the color mapping. Returns `(y, crange, s)`:
transformed values, range of the color scale (in transformed units) and scale `s`
(relevant only for `:asinh`). See the `plot` method of [`SecondaryStructureDiagram`](@ref)
for the meaning of the keywords.
"""
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
    _assign_secondary_structure(atoms, residues)

Compute the secondary structure with DSSP (`dssp_run`) and return one code per residue,
in the order of `residues`. If DSSP returns a different number of residues, match them
by (chain, residue number); residues without assignment are marked as loops (`"C"`).
"""
function _assign_secondary_structure(atoms, residues)
    ssdata = dssp_run(atoms)

    nres = length(residues)
    length(ssdata) == nres && return [string(s.sscode) for s in ssdata]

    s1 = first(ssdata)
    if hasproperty(s1, :chain) && hasproperty(s1, :resnum)
        lookup = Dict((string(s.chain), Int(s.resnum)) => string(s.sscode) for s in ssdata)
        codes = String[]
        n_unmatched = 0
        for r in residues
            a = first(r)
            code = get(lookup, (string(chain(a)), Int(resnum(a))), nothing)
            if isnothing(code)
                n_unmatched += 1
                code = "C"
            end
            push!(codes, code)
        end
        n_unmatched > 0 && @warn "$n_unmatched residue(s) without secondary structure assignment; marked as loop (C)."
        return codes
    end

    # COV_EXCL_START (defensive: dssp_run always returns chain and residue numbers)
    throw(ArgumentError("dssp_run returned $(length(ssdata)) residues, but the selection has $nres. " *
        "Check the selection (e.g. use \"protein\" to exclude water, ions and ligands)."))
    # COV_EXCL_STOP
end

"""
    secondary_structure_diagram(atoms::AbstractVector{<:Atom}, values::AbstractVector{<:Real})
    secondary_structure_diagram(atoms::AbstractVector{<:Atom}; property = beta, aggregate = mean)

Compute the secondary structure of each residue of `atoms` with DSSP (through
[`dssp_run`](@ref)) and associate one value to each residue. Returns a
[`SecondaryStructureDiagram`](@ref), which is drawn with `plot` (requires `CairoMakie`).

The values are given as a vector with one value per residue, in the order of
`eachresidue(atoms)`, or computed from a per-atom property:

- `property`: function applied to each atom (default: `beta`, the B-factor column).
- `aggregate`: how the values of the atoms of each residue are combined, called as
  `aggregate(property, residue)`. Default: the mean. Use `sum` for additive per-atom
  values (e.g. atomic contributions to a coordination number, from ComplexMixtures.jl),
  and the mean for intensive properties (B-factors), or when the same per-residue value
  was written in all atoms of the residue. `maximum` and `minimum` can also be used.

Helices shorter than 4 residues and strands shorter than 3 residues are converted to
loops.

# Example

```julia-repl
julia> using PDBTools, LazyArtifacts

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> diagram = secondary_structure_diagram(atoms)
```
"""
function secondary_structure_diagram(atoms::AbstractVector{<:Atom}, values::AbstractVector{<:Real})
    residues = collect(eachresidue(atoms))
    nres = length(residues)
    nres > 0 || throw(ArgumentError("The atom selection has no residues."))
    length(values) == nres ||
        throw(ArgumentError("values has $(length(values)) elements, but there are $nres residues."))

    chains = [string(chain(first(r))) for r in residues]
    residue_numbers = [Int(resnum(first(r))) for r in residues]
    segment_start = _segment_starts(chains, residue_numbers)

    # Secondary structure; short helices and strands become loops
    codes = _assign_secondary_structure(atoms, residues)
    for (i, j) in _structure_runs([c in _HELIX_CODES for c in codes], segment_start)
        (j - i + 1) < 4 && (codes[i:j] .= "C")
    end
    for (i, j) in _structure_runs([c in _STRAND_CODES for c in codes], segment_start)
        (j - i + 1) < 3 && (codes[i:j] .= "C")
    end

    labels = residue_ticks(atoms; oneletter = true)[2]
    length(labels) == nres || error("Residue labels ($(length(labels))) do not match residues ($nres).") # COV_EXCL_LINE

    return SecondaryStructureDiagram(string.(labels), chains, residue_numbers, codes, Float64.(values))
end

function secondary_structure_diagram(atoms::AbstractVector{<:Atom}; property = beta, aggregate = _residue_mean)
    residues = collect(eachresidue(atoms))
    # If the property is constant within each residue, it is probably already a
    # per-residue value; summing would multiply it by the number of atoms.
    if aggregate === sum
        per_atom = [property.(r) for r in residues]
        # Base.all: PDBTools defines its own `all` (the "all" selection)
        if Base.all(v -> Base.all(≈(first(v)), v), per_atom) && any(v -> length(v) > 1, per_atom)
            @warn "aggregate = sum, but the property is constant within each residue " *
                  "(it looks like a per-residue value). Consider the mean."
        end
    end
    values = [aggregate(property, r) for r in residues]
    return secondary_structure_diagram(atoms, values)
end

@testitem "secondary structure diagram helpers" begin
    using PDBTools
    runs = PDBTools._structure_runs
    @test runs(Bool[0, 1, 1, 0, 1], falses(5)) == [(2, 3), (5, 5)]
    # a new segment breaks a run
    @test runs(Bool[1, 1, 1, 1], Bool[1, 0, 1, 0]) == [(1, 2), (3, 4)]
    @test isempty(runs(falses(3), falses(3)))
    @test_throws ArgumentError runs(Bool[1, 0], Bool[1])

    @test PDBTools._segment_starts(["A", "A", "A", "B"], [1, 2, 4, 5]) == [true, false, true, true]

    @test PDBTools._sample_quantile([1, 2, 3, 4, 5], 0.5) == 3.0
    @test PDBTools._sample_quantile([1, 2, 3, 4], 0.5) == 2.5
    @test PDBTools._sample_quantile([7], 0.9) == 7.0
    @test_throws ArgumentError PDBTools._sample_quantile([1, 2], 1.5)

    cv = PDBTools._color_values
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

@testitem "secondary structure assignment" begin
    using PDBTools, LazyArtifacts
    pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    residues = collect(eachresidue(pdb))
    codes = PDBTools._assign_secondary_structure(pdb, residues)
    @test length(codes) == length(residues)
    # residues that are not in the structure given to DSSP are matched by
    # (chain, residue number); the missing one is marked as loop
    first10 = [at for r in residues[1:10] for at in r]
    codes = @test_logs (:warn, r"without secondary structure") match_mode = :any PDBTools._assign_secondary_structure(first10, residues[1:11])
    @test length(codes) == 11
    @test codes[end] == "C"
end

@testitem "secondary_structure_diagram" begin
    using PDBTools, LazyArtifacts
    atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    residues = collect(eachresidue(atoms))
    nres = length(residues)

    # values from the B-factor column: mean (default) and sum
    diagram = secondary_structure_diagram(atoms)
    @test diagram isa SecondaryStructureDiagram
    @test length(diagram.values) == length(diagram.secondary_structure) == length(diagram.chains) ==
          length(diagram.residue_numbers) == length(diagram.residue_labels) == nres
    @test diagram.values[1] ≈ sum(beta, residues[1]) / length(residues[1])
    @test secondary_structure_diagram(atoms; aggregate = sum).values[1] ≈ sum(beta, residues[1])
    @test secondary_structure_diagram(atoms; property = occup, aggregate = maximum).values[1] ≈
          maximum(occup, residues[1])
    # short helices and strands were converted to loops
    for (i, j) in PDBTools._structure_runs([c in PDBTools._HELIX_CODES for c in diagram.secondary_structure],
                                            PDBTools._segment_starts(diagram.chains, diagram.residue_numbers))
        @test j - i + 1 >= 4
    end

    # explicit values
    values = collect(range(-1.0, 1.0; length = nres))
    @test secondary_structure_diagram(atoms, values).values == values
    @test_throws ArgumentError secondary_structure_diagram(atoms, ones(nres - 1))
    @test_throws ArgumentError secondary_structure_diagram(atoms[1:0], Float64[])

    # warning when aggregate = sum but the property is constant within each residue
    atoms_const = deepcopy(atoms)
    for (i, r) in enumerate(eachresidue(atoms_const)), at in r
        at.beta = i
    end
    @test_logs (:warn, r"aggregate = sum") match_mode = :any secondary_structure_diagram(atoms_const; aggregate = sum)

    # show
    out = sprint(show, MIME"text/plain"(), diagram)
    @test occursin("SecondaryStructureDiagram with $nres residues", out)
end
