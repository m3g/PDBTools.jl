export ss

const _SS_HELIX_CODES = ("G", "H", "I", "P")
const _SS_BETA_CODES = ("E", "B")

"""
    _ss_runs(mask, newseg)

Return the contiguous segments `(i, j)` where `mask` is true. A segment is broken
when `newseg[k]` is true (start of a new chain or a gap in the numbering), so that
a helix or strand never crosses the boundary between two chains.
"""
function _ss_runs(mask::AbstractVector{Bool}, newseg::AbstractVector{Bool})
    length(mask) == length(newseg) || throw(ArgumentError("mask and newseg must have the same length"))
    runs = Tuple{Int,Int}[]
    n = length(mask)
    i = 1
    while i <= n
        if mask[i]
            j = i
            while j < n && mask[j+1] && !newseg[j+1]
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

# Sample quantile with linear interpolation (same definition as the default of Statistics.quantile)
function _ss_quantile(v::AbstractVector{<:Real}, p::Real)
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
    _ss_transform_contrib(x; mode = :robust, range_mode = :auto, scale = nothing,
                          limit = nothing, clip_quantile = 0.99, signed_tol = 1e-6)

Transform the per-residue values for the color mapping. Returns `(y, crange, s)`:
transformed values, range of the color scale (in transformed units) and scale `s`
(relevant only for `:asinh`).

Modes (`mode`):
- `:robust` (default): linear, with `lim` = quantile `clip_quantile` of `|x|`. Values
  beyond `lim` take the color of the extreme.
- `:linear`: linear, with `lim = max(|x|)`.
- `:asinh`: `asinh(x / s)`, linear for `|x| ≪ s` and logarithmic for `|x| ≫ s`; useful
  when the values span several orders of magnitude. `s = lim / 20` by default, or `scale`.

Range (`range_mode`):
- `:auto`: symmetric `(-lim, lim)` if some value is smaller than `-signed_tol * max(|x|)`;
  otherwise `(0, lim)`. The relative tolerance keeps numerical noise (e.g. -1e-12)
  from changing the whole scale.
- `:symmetric`: always `(-lim, lim)` (zero at the center of the color map).
- `:positive`: always `(0, lim)`; negative values take the minimum color.

`limit` and `scale` (in original units) fix `lim` and `s`. Use the same values in
all figures that will be compared.
"""
function _ss_transform_contrib(x::AbstractVector{<:Real};
    mode::Symbol = :robust,
    range_mode::Symbol = :auto,
    scale::Union{Nothing,Real} = nothing,
    limit::Union{Nothing,Real} = nothing,
    clip_quantile::Real = 0.99,
    signed_tol::Real = 1e-6,
)
    xf = float.(x)
    ref = maximum(abs, xf)

    # Limit in original units
    lim0 = if !isnothing(limit)
        float(limit)
    elseif mode == :robust
        q = _ss_quantile(abs.(xf), clip_quantile)
        iszero(q) ? ref : q # many zeros: use the maximum
    elseif mode in (:linear, :asinh)
        ref
    else
        throw(ArgumentError("normalize must be :robust, :linear or :asinh, got $mode"))
    end
    lim0 = (iszero(lim0) || !isfinite(lim0)) ? 1.0 : lim0 # all values are zero
    lim0 > 0 || throw(ArgumentError("limit must be positive, got $lim0"))

    # Type of range
    signed = if range_mode == :auto
        minimum(xf) < -signed_tol * ref
    elseif range_mode == :symmetric
        true
    elseif range_mode == :positive
        false
    else
        throw(ArgumentError("range_mode must be :auto, :symmetric or :positive, got $range_mode"))
    end

    # Transformation
    s = 1.0
    if mode == :asinh
        s = isnothing(scale) ? lim0 / 20 : float(scale)
        s > 0 || throw(ArgumentError("scale must be positive, got $s"))
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
    _ss_assign(atoms, residues)

Compute the secondary structure with DSSP (`dssp_run`) and return one code per residue,
in the order of `residues`. If DSSP returns a different number of residues, try to match
them by (chain, residue number); if that is not possible, throw an error instead of
silently filling the vector.
"""
function _ss_assign(atoms, residues)
    ssdata = dssp_run(atoms)

    nres = length(residues)
    length(ssdata) == nres && return [string(s.sscode) for s in ssdata]

    s1 = first(ssdata)
    if hasproperty(s1, :chain) && hasproperty(s1, :resnum)
        lookup = Dict((string(s.chain), Int(s.resnum)) => string(s.sscode) for s in ssdata)
        ssvec = String[]
        n_unmatched = 0
        for r in residues
            a = first(r)
            code = get(lookup, (string(chain(a)), Int(resnum(a))), nothing)
            if isnothing(code)
                n_unmatched += 1
                code = "C"
            end
            push!(ssvec, code)
        end
        n_unmatched > 0 && @warn "$n_unmatched residue(s) without secondary structure assignment; marked as loop (C)."
        return ssvec
    end

    throw(ArgumentError("dssp_run returned $(length(ssdata)) residues, but the selection has $nres. " *
        "Check the selection (e.g. use \"protein\" to exclude water, ions and ligands)."))
end

"""
    ss(atoms::AbstractVector{<:Atom}, contrib = nothing; kwargs...)
    ss(pdbfile::AbstractString, selection = "protein", contrib = nothing; kwargs...)

Draw a per-residue secondary structure diagram (loop/helix/β, computed with DSSP through
[`dssp_run`](@ref)) colored by a per-residue
property (for example, the coordination number from ComplexMixtures.jl) or, if
`contrib = nothing`, by the B-factor (beta) column, aggregated per residue according
to `beta_reduce`. Returns the `CairoMakie.Figure`.

!!! note
    This function requires `using CairoMakie`.

# Keyword arguments
- `beta_reduce` (used only when `contrib = nothing`):
  - `:mean` (default): mean over the atoms of the residue. Use when beta is a
    crystallographic B-factor or when the same per-residue value was written to all atoms.
  - `:sum`: sum over the atoms of the residue. Use when beta contains per-atom
    contributions (e.g. atomic coordination numbers from ComplexMixtures.jl), which are additive.
- `color`: color map. If `nothing`, uses `:balance` (diverging, zero = white) when the
  range is symmetric and `:tempo` (sequential) when it is `(0, lim)`.
- `normalize`: `:robust` (default; linear with limit at the quantile `clip_quantile` of
  |value|), `:linear` (limit = maximum) or `:asinh` (for values spanning several orders
  of magnitude).
- `range_mode`: `:auto` (default), `:symmetric` or `:positive`.
- `limit`, `scale`: limit of the color scale and scale of `:asinh`, in original units.
  Use the same values to compare several systems on the same scale.
- `clip_quantile = 0.99`: quantile used by `:robust`.
- `signed_tol = 1e-6`: relative tolerance to decide whether there are negative values (`:auto`).
- `per_chain = false`: divide the values of each chain by the maximum |value| of that chain.
- `chain_new_row = true`: each chain starts in a new row, identified on the y axis.
- `colorbar = true`, `colorbar_label = ""`: horizontal color bar below the rows.
- `output = nothing`: if given, the figure is saved to this path, always as SVG (vector);
  another extension (or none) is replaced by `.svg`.

Size and font are fixed so that all figures come out identical and ready for
publication: full page width (≈ 17 cm), 50 residues per row, 9 pt per residue,
8 pt Computer Modern font (if the CMU Serif font files are found). The height grows
with the number of rows.

# Example

```julia-repl
julia> using PDBTools, CairoMakie

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> fig = ss(atoms; colorbar_label = "B-factor", output = "ss.svg");
```
"""
function ss end

@testitem "ss diagram helpers" begin
    using PDBTools
    runs = PDBTools._ss_runs
    @test runs(Bool[0, 1, 1, 0, 1], falses(5)) == [(2, 3), (5, 5)]
    # a new segment breaks a run
    @test runs(Bool[1, 1, 1, 1], Bool[1, 0, 1, 0]) == [(1, 2), (3, 4)]
    @test isempty(runs(falses(3), falses(3)))

    @test PDBTools._ss_quantile([1, 2, 3, 4, 5], 0.5) == 3.0
    @test PDBTools._ss_quantile([1, 2, 3, 4], 0.5) == 2.5
    @test PDBTools._ss_quantile([7], 0.9) == 7.0

    tc = PDBTools._ss_transform_contrib
    x = collect(0.0:10.0)
    y, crange, s = tc(x; mode = :linear)
    @test y == x
    @test crange == (0.0, 10.0)
    # negative values give a symmetric range; numerical noise does not
    @test tc([-5.0, 10.0]; mode = :linear)[2] == (-10.0, 10.0)
    @test tc([-1e-12, 10.0]; mode = :linear)[2] == (0.0, 10.0)
    @test tc(x; mode = :linear, range_mode = :symmetric)[2] == (-10.0, 10.0)
    # robust limit is the quantile of |x|
    @test tc(x; mode = :robust, clip_quantile = 0.5)[2] == (0.0, 5.0)
    # fixed limit; all-zero data do not break the scale
    @test tc(x; limit = 2.0)[2] == (0.0, 2.0)
    @test tc(zeros(4))[2] == (0.0, 1.0)
    # asinh
    y, crange, s = tc(x; mode = :asinh, scale = 1.0)
    @test y ≈ asinh.(x)
    @test crange[2] ≈ asinh(10.0)
    @test s == 1.0
    @test_throws ArgumentError tc(x; mode = :log)
    @test_throws ArgumentError tc(x; range_mode = :negative)
end
