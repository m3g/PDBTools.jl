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
SASA, RMSF, interaction energies, conservation scores, etc. The resulting
[`HelicalWheel`](@ref) is plotted with `plot` (a Makie backend must be loaded, e.g.
`using CairoMakie`):

```@example wheel
using PDBTools, CairoMakie
sequence = "LKKLLKLLKKLLKL"
values = [1.0, 0.1, 0.2, 0.9, 0.8, 0.1, 0.7, 1.0, 0.2, 0.1, 0.9, 0.8, 0.2, 0.9]
wheel = helical_wheel(sequence, values)
```

```@example wheel
plot(wheel; colorbar_title = "property")
```

The `HelicalWheel` contains the mean hydrophobicity of the helix and its hydrophobic
moment, in the `mean_hydrophobicity` and `μ` fields:

```@example wheel
wheel.μ
```

Some details of the drawing:

- the residues are placed at 100° from each other; after 18 residues (5 turns) a new
  ring starts, with a larger radius, so that long helices do not overlap;
- the arrow is the hydrophobic moment, and the wheel is rotated so that it always
  points down: the hydrophobic face of the helix is at the bottom of the figure;
- the first and last residues are labeled with their residue numbers (keyword
  `residue_numbers` of `helical_wheel`); use `residue_numbers = :all` (or `:none`) in
  `plot` to label all residues (or none);
- by default the values are normalized to `[0, 1]` (`normalize = true`). To compare
  several helices on the same color scale, use `normalize = false` and fix the limits
  with `clims`.

## Per-residue values

The values of the residues are the input of the helical wheel, and can come from any
source, as long as there is one value per residue, in the order of the sequence.

Here, we compute the solvent accessible surface area (SASA) of the atoms of a protein,
and the SASA of each residue of an α-helix of the protein (residues 13 to 29):

```@example wheel
using LazyArtifacts
pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
atomic_sasa = sasa_particles(pdb)
helix = select(pdb, "resnum >= 13 and resnum <= 29")
values = [ sasa(atomic_sasa, at -> residue(at) == residue(r)) for r in eachresidue(helix) ]
```

```@example wheel
sequence = join(oneletter.(resname.(eachresidue(helix))))
plot(helical_wheel(sequence, values); normalize = false, colorbar_title = "SASA (Å²)")
```

The buried face of the helix, which packs against the rest of the protein, is at the bottom.

Values computed by other packages can be used as well. For example, with a
`ResidueContributions` object `rc` from [ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl),
the contribution of each residue to the solvation of the protein at a distance `d` is:

```julia
id = findfirst(>=(d), rc.d)
values = [rc.residue_contributions[i][id] for i in eachindex(rc.residue_contributions)]
plot(helical_wheel(sequence, values); colorbar_title = "CN")
```

## From a structure

With a selection of atoms (of a single chain), the sequence, the residue numbers and the
values are taken from the structure itself. The values are computed from a per-atom
property (`property`, by default the `beta` field of the atoms), combined for each
residue by `aggregate`, as `aggregate(property, residue)`:

- use `aggregate = sum` for additive properties, such as the SASA, or atomic contributions
  to a coordination number;
- use the mean (default) for intensive properties, such as B-factors or RMSF.

Residues (and chains, segments, models and structures) are vectors of atoms, such that
`aggregate(property, residue)` is, for example, `sum(beta, residue)` or `mean(beta, residue)`
(with `using Statistics`).

Per-atom values computed by other programs are often stored in the `beta` column of PDB files.
Here, we store the SASA of each atom in the `beta` field of the atoms, and the helical wheel
above is obtained directly from the helix:

```@example wheel
for (i, atom) in enumerate(pdb)
    atom.beta = atomic_sasa[i] # SASA of atom i (the atoms of helix are those of pdb)
end
wheel = helical_wheel(helix; aggregate = sum)
plot(wheel; normalize = false, colorbar_title = "SASA (Å²)")
```

The `property` keyword selects another per-atom property, for example
`property = occup` for the occupancy column. The first and last residues are labeled
with their residue numbers in the structure.

## Hydrophobic moment and mean hydrophobicity

[`hydrophobic_moment`](@ref) computes the hydrophobic moment (Eisenberg) of a sequence,
without drawing, as a 2D vector, whose magnitude is the amphipathicity of the helix:

```@example wheel
m = hydrophobic_moment("LKKLLKLLKKLLKL")
```

```@example wheel
hypot(m...)
```

[`mean_hydrophobicity`](@ref) computes the mean hydrophobicity of the residues of a sequence:

```@example wheel
mean_hydrophobicity("LKKLLKLLKKLLKL")
```

The default hydrophobicity scale is that of Fauchère & Pliska (1983)[^1]. Another scale can
be given as a dictionary from one-letter codes to values, with the `hydrophobicity` keyword
of `hydrophobic_moment`, `mean_hydrophobicity` and `helical_wheel`.

[^1]: J. Fauchère and V. Pliska, *Hydrophobic parameters π of amino-acid side chains from the partitioning of N-acetyl-amino-acid amides*, Eur. J. Med. Chem. **18**, 369-375 (1983).

## Reference

```@docs
helical_wheel
HelicalWheel
plot(::HelicalWheel)
hydrophobic_moment
mean_hydrophobicity
```
