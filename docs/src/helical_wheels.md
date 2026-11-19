```@meta
CollapsedDocStrings = true
```

# Helical wheels

A helical wheel is a projection of an α-helix along its axis. Consecutive residues
are separated by 100° (3.6 residues per turn), so residues that face the same side of
the helix appear next to each other in the wheel. Coloring the residues by a
per-residue property shows at a glance which face of the helix carries that property,
for example, which face interacts with the membrane, the solvent or another protein.

## Drawing a helical wheel

[`helical_wheel`](@ref) receives the sequence of the helix (one-letter codes) and one
value per residue, in the same order. Any property can be used: coordination numbers,
SASA, RMSF, interaction energies, conservation scores, etc. `Plots` must be loaded.

```@example wheel
using PDBTools, Plots
sequence = "LKKLLKLLKKLLKL"
values = [1.0, 0.1, 0.2, 0.9, 0.8, 0.1, 0.7, 1.0, 0.2, 0.1, 0.9, 0.8, 0.2, 0.9]
p, info = helical_wheel(sequence, values; colorbar_title = "property")
p
```

The function returns the plot and a named tuple with the mean hydrophobicity of the
helix and its hydrophobic moment:

```@example wheel
info
```

Some details of the drawing:

- the residues are placed at 100° from each other; after 18 residues (5 turns) a new
  ring starts, with a larger radius, so that long helices do not overlap;
- the arrow is the hydrophobic moment, and the wheel is rotated so that it always
  points down: the hydrophobic face of the helix is at the bottom of the figure;
- the first and last residues are labeled with their residue numbers
  (keyword `residue_numbers`);
- by default the values are normalized to `[0, 1]` (`normalize = true`). To compare
  several helices on the same color scale, use `normalize = false` and fix the limits
  with `clims`.

## From a structure

With a selection of atoms, the sequence, the residue numbers and the values are taken
from the structure itself. The values are computed from a per-atom property (by
default, the `beta` column of the PDB file) combined for each residue by `aggregate`:

- use `aggregate = sum` for additive properties, such as atomic contributions to a
  coordination number written in the `beta` column;
- use the mean (default) for intensive properties, such as B-factors or RMSF.

The selection must contain a single chain.

```@example wheel
pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
first_atom = first(pdb)
helix = select(pdb, at -> chain(at) == chain(first_atom) && resnum(at) < resnum(first_atom) + 18)
p, info = helical_wheel(helix; colorbar_title = "B-factor")
p
```

The `property` keyword selects another per-atom property, for example
`property = occup` for the occupancy column.

## Values from other packages

The values can come from any source. For example, with a `ResidueContributions` object
`rc` from [ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl), the
contribution of each residue at a distance `d` is:

```julia
id = findfirst(>=(d), rc.d)
values = [rc.residue_contributions[i][id] for i in eachindex(rc.residue_contributions)]
p, info = helical_wheel(sequence, values; colorbar_title = "CN")
```

## Hydrophobic moment and per-residue values

[`hydrophobic_moment`](@ref) computes the mean hydrophobicity and the hydrophobic
moment (Eisenberg) of a sequence, without drawing. The default hydrophobicity scale
is that of Fauchère & Pliska; another scale can be given as a dictionary from
one-letter codes to values, with the `hydrophobicity` keyword.

```@example wheel
hydrophobic_moment("LKKLLKLLKKLLKL")
```

[`residue_values`](@ref) returns the per-residue values, the sequence and the residue
numbers of a selection, as used by `helical_wheel`:

```@example wheel
values, sequence, numbers = residue_values(helix; aggregate = sum)
sequence
```

## Reference

```@docs
helical_wheel
hydrophobic_moment
residue_values
PDBTools.FAUCHERE_PLISKA
```
