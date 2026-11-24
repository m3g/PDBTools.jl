```@meta
CollapsedDocStrings = true
```

# Secondary structure diagrams

A secondary structure diagram shows the secondary structure of each residue along the
sequence, with each residue colored by a value: helices as ribbons, β strands as arrows
and loops as thin bars. [`secondary_structure_diagram`](@ref) computes the secondary structure and the
values of the residues, and the resulting [`SecondaryStructureDiagram`](@ref) is plotted with `plot`
(a Makie backend must be loaded, e.g. `using CairoMakie`).

## Coloring by a property of the atoms

By default, the residues are colored by the mean B-factor of their atoms (in this test
structure, the B-factors grow from the core to the surface of the protein):

```@example ss-diagram
using PDBTools, LazyArtifacts, CairoMakie
atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
d = secondary_structure_diagram(atoms)
```

```@example ss-diagram
plot(d; colorbar_title = "B-factor")
```

The secondary structure is computed with DSSP (`ss_method = :dssp`, default) or STRIDE
(`ss_method = :stride`). Helices shorter than 4 residues and strands shorter than 3
residues are drawn as loops. Each row has 50 residues, and each chain starts a new row,
labeled on the left (`chain_new_row = true`).

The value of each residue is computed from a per-atom property (`property`, by default the
`beta` field of the atoms), combined for each residue by `aggregate`, as
`aggregate(property, residue)`:

- the mean (default): for intensive properties, such as B-factors or RMSF;
- `aggregate = sum`: for additive per-atom values, such as atomic contributions to a
  coordination number from [ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl).

`property` can be any function of the atom. Here, each residue is colored by the mean
distance of its atoms to the center of mass:

```@example ss-diagram
c = center_of_mass(atoms)
distance_to_center(at) = sqrt((at.x - c[1])^2 + (at.y - c[2])^2 + (at.z - c[3])^2)
plot(secondary_structure_diagram(atoms; property = distance_to_center); colorbar_title = "distance to center (Å)")
```

## Coloring by custom values

The values can also be given directly, one per residue, in the order of
`eachresidue(atoms)`:

```@example ss-diagram
nres = length(eachresidue(atoms))
values = [sin(i / 10) for i in 1:nres]
plot(secondary_structure_diagram(atoms, values); colorbar_title = "sin of residue index")
```

## Using a precomputed secondary structure

To avoid computing the secondary structure again, the output of [`dssp_run`](@ref) or
[`stride_run`](@ref) can be given directly. Without values, only the secondary structure
is drawn:

```@example ss-diagram
ss = dssp_run(atoms)
plot(secondary_structure_diagram(ss))
```

The values are then given as a second argument, one per element of the output of 
`dssp_run` or `stride_run`, which has one element per protein residue (residues that DSSP
or STRIDE cannot classify, e.g. with missing backbone atoms, are assigned as loop or coil):

```@example ss-diagram
values = [sum(beta, r) / length(r) for r in eachresidue(atoms)]
plot(secondary_structure_diagram(ss, values); colorbar_title = "B-factor")
```

## Saving

The figure has a fixed layout (rows of 50 residues, ≈ 24 cm wide, and a 12 pt font), so
that all diagrams come out identical. To keep the sizes in typographic points, save it as
SVG with CairoMakie and `pt_per_unit = 1` (with another backend, call 
`CairoMakie.activate!()` before saving). With `pt_per_unit = 2/3`, the figure has the width
of a full page (≈ 16 cm), with an 8 pt font:

```julia
fig = plot(d; colorbar_title = "B-factor")
CairoMakie.save("diagram.svg", fig; pt_per_unit = 1) # `save` is also exported by PDBTools
CairoMakie.save("diagram_page_width.svg", fig; pt_per_unit = 2/3)
```

The background is white by default; use `background_color = :transparent` for a
transparent background.

To set the size of the figure, use `size = (width, height)` (in pt). The rows are then
stretched to fill the figure, and the font size is kept:

```julia
fig = plot(d; size = (800, 600))
```

## Color scale

When some values are negative, the color scale is symmetric around zero and uses a
diverging color map (`:balance`); otherwise it goes from zero to the maximum and uses a
sequential color map (`:tempo`). This can be changed with `range_mode`
(`:auto`, `:symmetric` or `:positive`) and `colormap`.

The keyword `color_scale` defines how the values are mapped to colors:

- `:robust` (default): linear, with the upper limit at the 99% quantile of the absolute
  values (`clip_quantile`), so that a few extreme residues do not wash out the others;
- `:linear`: linear, with the upper limit at the maximum absolute value;
- `:asinh`: linear for small values and logarithmic for large ones, for values that
  span several orders of magnitude.

To compare several structures on the same scale, fix the limit (and, for `:asinh`, the
scale) in the units of the values, with the same values in all figures:

```julia
fig1 = plot(secondary_structure_diagram(atoms1, values1); color_limit = 2.0)
fig2 = plot(secondary_structure_diagram(atoms2, values2); color_limit = 2.0)
```

With `normalize_per_chain = true`, the values of each chain are divided by the maximum
absolute value of that chain.

## Reference

```@docs
secondary_structure_diagram
SecondaryStructureDiagram
plot(::SecondaryStructureDiagram)
```
