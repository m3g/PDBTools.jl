export plot_ss_diagram

"""
    plot_ss_diagram(atoms::AbstractVector{<:Atom}; property = :beta, aggregate = mean, ss_method = :dssp, kwargs...)
    plot_ss_diagram(ss::AbstractVector{<:SSData}; property = nothing, kwargs...)

Draw the secondary structure of each residue along the sequence (helices as ribbons,
β strands as arrows and loops as thin bars), with each residue colored by a value.
Requires `using CairoMakie`, and returns a `CairoMakie.Figure`.

# Input

- `atoms`: the secondary structure is computed with `dssp_run` (`ss_method = :dssp`,
  default) or `stride_run` (`ss_method = :stride`).
- `ss`: the output of [`dssp_run`](@ref) or [`stride_run`](@ref), a `Vector{SSData}`,
  so that the secondary structure is not computed again.

Helices shorter than 4 residues and strands shorter than 3 residues are drawn as loops.

# Values of the residues

With `atoms`, `property` can be:
- a `Symbol` or `String` with the name of an atom field, combined for the atoms of each
  residue by `aggregate`: `:beta` (default, the B-factor column), `:occupancy` (or `:occup`),
  `:charge`, etc.;
- a function of the atom, e.g. `at -> at.beta^2`, also combined by `aggregate`;
- a vector with one value per residue, in the order of `eachresidue(atoms)`.

`aggregate` is called as `aggregate(f, residue)`. Default: the mean. Use `sum` for additive
per-atom values (e.g. atomic contributions to a coordination number, from ComplexMixtures.jl).
`maximum` and `minimum` can also be used.

With `ss`, `property` is `nothing` (default: only the drawing, without colors) or a vector
with one value per element of `ss`.

# Layout

The figure has a fixed layout, so that all diagrams come out identical: 50 residues per
row, 13.5 pt per residue (≈ 24 cm wide rows) and a 12 pt font (Computer Modern, `cmr10`,
if available). The height grows with the number of rows. To keep the sizes in typographic
points, save as SVG with `CairoMakie.save("diagram.svg", fig; pt_per_unit = 1)` (`save`
alone is ambiguous, because it is also exported by PDBTools). With `pt_per_unit = 2/3`,
the figure has the width of a full page (≈ 16 cm), with an 8 pt font.

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
- `font = "cmr10"`: font of the labels. Makie uses a default font if it is not found.
- `size = nothing`: size of the figure, as `(width, height)`, in pt. By default, the
  size is fixed by the content (675 pt wide rows, 60 pt high), such that all diagrams
  have the same scale. With a given size, the rows are stretched to fill the figure
  (the font size is not changed).

Values beyond the color limits take the color of the extreme of the color map.

# Examples

```julia-repl
julia> using PDBTools, LazyArtifacts, CairoMakie

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> fig = plot_ss_diagram(atoms; property = :beta, colorbar_title = "B-factor");

julia> ss = dssp_run(atoms);

julia> fig = plot_ss_diagram(ss); # only the drawing

julia> CairoMakie.save("diagram.svg", fig; pt_per_unit = 1)
```
"""
function plot_ss_diagram end
