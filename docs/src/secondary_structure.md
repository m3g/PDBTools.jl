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

[`secondary_structure_diagram`](@ref) computes the secondary structure of each residue
with DSSP (through [`dssp_run`](@ref)) and associates one value to each residue. The
result, a [`SecondaryStructureDiagram`](@ref), is drawn with `plot`, which shows the
secondary structure along the sequence with each residue colored by its value.
`CairoMakie` must be loaded.

```@example ss-diagram
using PDBTools, LazyArtifacts, CairoMakie
atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
diagram = secondary_structure_diagram(atoms)
```

```@example ss-diagram
plot(diagram; colorbar_title = "B-factor")
```

Helices are drawn as ribbons, β strands as arrows and loops as thin bars. Helices
shorter than 4 residues and strands shorter than 3 residues are drawn as loops. Each
row has 50 residues, and each chain starts a new row, labeled on the left
(`chain_new_row = true`).

The figure has a fixed size (full page width, about 17 cm) and a fixed font size
(8 pt), so that all diagrams come out identical and ready for publication. If the
CMU Serif (Computer Modern) font files are installed, they are used. To keep the sizes
in typographic points, save the figure as SVG with `pt_per_unit = 1`:

```julia
fig = plot(diagram; colorbar_title = "B-factor")
save("diagram.svg", fig; pt_per_unit = 1)
```

If `Plots` is also loaded, call `CairoMakie.plot(diagram)`.

### Values of each residue

Without a second argument, the value of each residue is computed from a per-atom
property, by default the `beta` column of the structure, combined for each residue by
`aggregate`:

- the mean (default): for B-factors, or when the same per-residue value was written in
  all atoms of the residue;
- `aggregate = sum`: for additive per-atom values, such as atomic contributions to a
  coordination number from [ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl).

Another per-atom property can be chosen with `property`, for example
`property = occup`. Any per-residue values can also be given as a vector, with one value
per residue, in the order of `eachresidue(atoms)`:

```@example ss-diagram
nres = length(eachresidue(atoms))
values = [sin(i / 10) for i in 1:nres]
plot(secondary_structure_diagram(atoms, values); colorbar_title = "property")
```

### Color scale

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
fig1 = plot(diagram1; color_limit = 2.0)
fig2 = plot(diagram2; color_limit = 2.0)
```

With `normalize_per_chain = true`, the values of each chain are divided by the maximum
absolute value of that chain.

```@docs
secondary_structure_diagram
SecondaryStructureDiagram
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




