```@meta
CollapsedDocStrings = true
```

# [Secondary structure](@id secondary-structure)

These functions provide an interface to compute the secondary structure assignment of proteins
using the STRIDE and DSSP algorithms, using as inputs vectors of `PDBTools.Atom`s.

```@docs
stride_run
dssp_run
```

The `stride_run` and `dssp_run` functions return a vector of `SSData` objects, each containing the secondary structure assignment and backbone dihedral angles for each residue.

These functions return a vector of `SSData` objects, as defined in [ProteinSecondaryStructures.jl](https://biojulia.dev/ProteinSecondaryStructures.jl). The secondary structure assignment codes are available 
in the [ProteinSecondaryStructures.jl documentation](https://biojulia.dev/ProteinSecondaryStructures.jl/stable/overview/#Secondary-structure-classes).

!!! note
    Non-standard residue or atom names may lead to incorrect secondary structure assignments. To ensure accurate results, it is recommended to replace non-standard names with their standard three-letter codes before running these functions, and to verify that all backbone atoms (N, CA, C, O) are present in each residue. Errors or warnings will be issued if residue names are not recognized or backbone atoms are missing.

## Example using STRIDE

The `stride_run` function runs the STRIDE algorithm on the provided array of atoms.
```@example secondary-structure
using PDBTools, LazyArtifacts
atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
ss = stride_run(atoms)
```

To run with DSSP just use the `dssp_run` function instead. 

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

## Utility functions 

`PDBTools` also reexports the `ss_composition` , `ss_name`, `ss_code`, `ss_number` of `ProteinSecondaryStructures.jl`, which can be used to analyze the secondary structure assignment results. For example, to compute the secondary structure composition from the `ss` vector obtained above:

```@example secondary-structure
ss_composition(ss)
```

```@example secondary-structure
ss_name(ss[1]) # name of the secondary structure of the first residue
```

For further information, refer to the [ProteinSecondaryStructures.jl documentation](https://biojulia.dev/ProteinSecondaryStructures.jl/stable/user_guide/#Retrieving-names,-codes,-and-numeric-codes).




