```@meta
CollapsedDocStrings = true
```
# Dihedrals and Ramachandran plots

Dihedral angles can be computed with the `dihedral` function, and an application 
of this function is the computation of Ramachandran plots. 

```@docs
dihedral
Ramachandran
scatter(::Ramachandran)
```

## Dihedral angles

The `dihedral` function computes the dihedral angle given four atoms:

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
prot = read_pdb(PDBTools.TESTPDB, "protein");
ram = Ramachandran(prot)
```

Given the `ram::Ramachandran` object, the `scatter` function from `Plots` can be used to 
produce the Ramachandran plot:

```@example ramachandran
using Plots
scatter(ram)
```

All `scatter` parameters can be customized using the `Plots` keyword syntax. 

## Setting backbone dihedral angles

The `set_phi!` and `set_psi!` functions set the phi and psi backbone dihedral angles of a
protein to given target values, rigidly rotating the corresponding downstream fragment of
the chain (side chain, and every atom from that point to the C-terminus) around the
relevant backbone bond. Rotations are applied residue by residue, from the first to the
last, so setting phi does not disturb any previously-set psi angle, and vice versa.

```@docs
set_phi!
set_psi!
```

Both functions require `atoms` to belong to a single protein chain, with the backbone atoms
needed to define every dihedral to be set. Angles are given in radians by default:

```@example ramachandran
using PDBTools
prot = read_pdb(PDBTools.TESTPDB, "protein");
nres = length(collect(eachresidue(prot)));
set_phi!(prot, fill(-60.0, nres - 1); unit="deg")
set_psi!(prot, fill(-45.0, nres - 1); unit="deg")
ram = Ramachandran(prot)
```

## Fully extended chain

`extended_chain!` sets every phi and psi angle of a protein to 180°, producing the fully
extended (all-trans, "C5") backbone conformation - for example, to use as an explicit
model of the unfolded state (as required by the [`MTRecord`](@ref record_model) transfer
model). `extended_chain` is the non-mutating version, returning a new, independent vector
of atoms:

```@docs
extended_chain!
extended_chain
```

```@example ramachandran
ext = extended_chain(prot)
ram_ext = Ramachandran(ext)
```

## Check the stereochemistry of protein residues

```@docs
zeta
zeta_check
```
