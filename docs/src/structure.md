# Structures and unit cells

`read_pdb`, `read_mmcif` and `wget` return a `Structure`: a vector of atoms that also carries
the metadata of the system, in particular its unit cell. A `Structure` behaves as a vector of 
`Atom`s, such that all functions that operate on vectors of atoms also operate on structures.

```@docs
Structure
PDBTools.MetaData
```

## Reading structures

The unit cell is read from the `CRYST1` field of PDB files, or from the `_cell` fields of mmCIF
files:

```@example structure
using PDBTools
str = read_pdb(PDBTools.test_file("pbc.pdb"))
```

```@example structure
str.unitcell
```

If the file does not define a unit cell (or defines the 1Å unit cube, which indicates that the
structure is not crystallographic, as in NMR or cryo-EM structures), `str.unitcell` is `nothing`.

## Indexing and iteration

A `Structure` can be indexed and iterated as a vector of atoms:

```@example structure
str[1]
```

```@example structure
collect(eachchain(str))
```

Indexing with ranges, selections, `filter`, `copy`, and other functions that create new vectors of 
atoms return new `Structure`s, with a copy of the metadata. Thus, the unit cell is preserved:

```@example structure
protein = select(str, "protein")
protein.unitcell
```

The underlying vector of atoms is available as `str.atoms`.

## Unit cells

The unit cell is set with a [`UnitCell`](@ref), which can be constructed from the lattice 
parameters, from the sides of an orthorhombic box, or from a matrix with the lattice vectors 
as columns:

```@example structure
str.unitcell = UnitCell(85.0, 85.0, 85.0)
```

```@example structure
str.unitcell = UnitCell(a=85.0, b=85.0, c=85.0, γ=120.0)
```

```@docs
UnitCell
```

Unit cells are written by `write_pdb` and `write_mmcif`. 

## Periodic boundary conditions

The functions that support periodic boundary conditions (for example, [`contact_map`](@ref), 
[`hydrogen_bonds`](@ref), [`sasa_particles`](@ref), and [`visualize`](@ref)) accept the 
`unitcell` keyword argument as:

- `unitcell=false` (default): no periodic boundary conditions.
- `unitcell=true`: use the unit cell of the `Structure`.
- `unitcell=UnitCell(...)`: use the given unit cell.

```@example structure
str = read_pdb(PDBTools.test_file("pbc.pdb"), "protein")
sasa(sasa_particles(str; unitcell=true))
```
