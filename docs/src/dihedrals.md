```@meta
CollapsedDocStrings = true
```
# Dihedrals and Ramachandran plots

Dihderal angles can be computed with the `dihedral` function, and an application 
of this function is the computation of Ramachandran plots. 

```@docs
dihedral
Ramachandran
scatter(::Ramachandran)
```

## Dihedral angles

The dihderal function computes the dihderal angle given four atoms:

```jldoctest
julia> using PDBTools

julia> prot = read_pdb(PDBTools.TESTPDB, "protein");

julia> dihedral(prot[1], prot[5], prot[11], prot[13])
64.07296f0
```

## Ramachandran plot

The `Ramachandran` function and object are used to compute and plot Ramachandran plots for a protein structure. The call to `Ramachandran(vec)` where `vec` is a vector of `Atom`s returns a `Ramachandran` object, with fields `phi` and `psi`, containing the list of dihedral angles:

```@example ramachandran
using PDBTools
prot = read_pdb(PDBTools.TESTPDB);
ram = Ramachandran(prot)
```

Given the `ram::Ramachandran` object, the `scatter` function from `Plots` can be used to 
produce the Ramachandran plot:

```@example ramachandran
using Plots
scatter(ram)
```

All `scatter` parametes can be customized using the `Plots` keyword syntax. 