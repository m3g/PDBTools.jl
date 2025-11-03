```@meta
CollapsedDocStrings = true
```

# Solvent Accessible Surface Area (SASA)

These functions are used to compute the solvent accessible surface area (SASA) of structures or parts of a structure. They provide a very fast implementation of the [Shake-Rupley](https://doi.org/10.1016/0022-2836(73)90011-9) method, using a Fibbonacy lattice to construct the grid points.

```@docs
sasa_particles
sasa
```

A typical run of these functions consists in providing the structure of a protein to the first function, `sasa_particles`, to obtain a `SASA` object, which contains the accessible area per atom:

```@example sasa
using PDBTools
prot = read_pdb(PDBTools.TESTPDB, "protein")
atom_sasa = sasa_particles(prot)
```

The `atom_sasa` object created above can be used to extract the total accessible area or the accessible area of any sub-surface. The `sasa` function provides an interface for those extractions:

```@example sasa
sasa(atom_sasa) # total
```
```@example sasa
sasa(atom_sasa, "polar") 
```

```@example sasa
sasa(atom_sasa, "backbone")
```

```@example sasa
sasa(atom_sasa, "resname THR and residue < 50") 
```

In some situations, it might be useful to visualize the surface. The dots that form the surface can be obtained by running `sasa_particles` with the `output_dots` option set to `true`. Here, we use fewer dots for better visualization:

```@example sasa
atom_sasa = sasa_particles(prot; n_dots=100, output_dots=true) 
```

Where the `atom_sasa.dots` field contais the dots that are accessible to the surface for each atom. These can be plotted, for example, with:
```@example sasa
using Plots
dots = reduce(vcat, atom_sasa.dots)
scatter(Tuple.(coor.(prot)); color=:orange, msw=0, label="") # atom coordinates
scatter!(Tuple.(dots); # surface dots
    color=:blue, ms=1, msw=0, ma=0.5, # marker properties
    label="",
)
```

!!! note
    The `sasa_particles` function supports periodic boundary conditions if a unit cell is provided. 
    See the how to [read the unitcell](@ref read-unitcell)  for further information.
