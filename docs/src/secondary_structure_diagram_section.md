## Secondary structure diagrams

[`ss`](@ref) draws the secondary structure of each residue (loop, helix or β strand,
computed with DSSP through [`dssp_run`](@ref)) as a diagram along the sequence, with each
residue colored by a per-residue property. `CairoMakie` must be loaded.

```@example ss-diagram
using PDBTools, CairoMakie
atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
ss(atoms; colorbar_label = "B-factor")
```

Helices are drawn as a ribbon, β strands as arrows and loops as thin bars. Helices
shorter than 4 residues and strands shorter than 3 residues are drawn as loops. Each
row has 50 residues, and each chain starts a new row, labeled on the left
(`chain_new_row = true`).

The figure has a fixed size (full page width, about 17 cm) and a fixed font size
(8 pt), so that all diagrams come out identical and ready for publication. If the
CMU Serif (Computer Modern) font files are installed, they are used. With the
`output` keyword, the figure is saved as SVG (vector):

```julia
ss(atoms; colorbar_label = "B-factor", output = "ss.svg")
```

### Coloring by a per-residue property

Without a second argument, the residues are colored by the `beta` column of the
structure, combined for each residue by `beta_reduce`:

- `:mean` (default): for B-factors, or when the same per-residue value was written in
  all atoms of the residue;
- `:sum`: for additive per-atom values, such as atomic contributions to a coordination
  number from [ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl).

Any other per-residue property can be given as a vector, with one value per residue
of the selection:

```@example ss-diagram
nres = length(eachresidue(atoms))
values = [sin(i / 10) for i in 1:nres]
ss(atoms, values; colorbar_label = "property")
```

The structure can also be read directly from a file (PDB or mmCIF), with a selection:

```julia
ss("structure.pdb", "protein and chain A", values; output = "ss.svg")
```

### Color scale

When some values are negative, the color scale is symmetric around zero and uses a
diverging color map (`:balance`); otherwise it goes from zero to the maximum and uses a
sequential color map (`:tempo`). This can be changed with `range_mode`
(`:auto`, `:symmetric` or `:positive`) and `color`.

The keyword `normalize` defines how the values are mapped to colors:

- `:robust` (default): linear, with the upper limit at the 99% quantile of the absolute
  values (`clip_quantile`), so that a few extreme residues do not wash out the others;
- `:linear`: linear, with the upper limit at the maximum absolute value;
- `:asinh`: linear for small values and logarithmic for large ones, for values that
  span several orders of magnitude.

To compare several structures on the same scale, fix the limit (and, for `:asinh`, the
scale) in original units, with the same values in all figures:

```julia
ss(atoms1, values1; limit = 2.0, output = "system1.svg")
ss(atoms2, values2; limit = 2.0, output = "system2.svg")
```

With `per_chain = true`, the values of each chain are divided by the maximum absolute
value of that chain.

```@docs
ss
```
