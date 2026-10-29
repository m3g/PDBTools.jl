# Structures and unit cells

`read_pdb`, `read_mmcif`, `read_xyz` and `wget` return a `Structure`: a vector of atoms that 
also carries the metadata of the system: its unit cell, title, and connectivity. A `Structure` 
behaves as a vector of `Atom`s, such that all functions that operate on vectors of atoms also 
operate on structures.

```@docs
Structure
PDBTools.MetaData
```

## [Metadata](@id structure-metadata)

The metadata fields are accessible as properties of the structure, and are `nothing` when not 
defined:

| Property | Read from | Written to |
|:---------|:----------|:-----------|
| `unitcell` | `CRYST1` (PDB), `_cell` (mmCIF) | `CRYST1` (PDB), `_cell` (mmCIF) |
| `title` | `TITLE` (PDB), `_struct.title` (mmCIF), comment line (xyz) | the same fields |
| `connectivity` | `CONECT` (PDB), `_geom_bond` (mmCIF) | `CONECT` (PDB), `_geom_bond` (mmCIF) |

The metadata can be set in the `Structure` constructor (e.g. `Structure(atoms; title="My system")`),
or by assigning the properties (e.g. `str.title = "My system"`).

## Reading structures

The unit cell is read from the `CRYST1` field of PDB files, or from the `_cell` fields of mmCIF
files:

```@example structure
using PDBTools, LazyArtifacts
str = read_pdb(PDBTools.test_file("pbc.pdb"))
```

```@example structure
str.unitcell
```

If the file does not define a unit cell (or defines the 1Å unit cube, which indicates that the
structure is not crystallographic, as in NMR or cryo-EM structures), `str.unitcell` is `nothing`.

## Title

The title is read from the `TITLE` records of PDB files, the `_struct.title` field of mmCIF files, 
or the comment line of xyz files:

```@example structure
ethanol = read_pdb(PDBTools.test_file("ethanol.pdb"))
ethanol.title
```

## Connectivity

The bonds of the atoms are read from the `CONECT` records of PDB files, or from the `_geom_bond`
fields of mmCIF files. The connectivity is a `Vector{Vector{Int32}}`, where `connectivity[i]` 
contains the indices of the atoms bonded to the atom with index `i`:

```@example structure
ethanol.connectivity
```

The indices are the sequential indices of the atoms in the file, `index(atom)`, which are not 
changed by selections, filtering, or reordering of the atoms. Thus, the connectivity remains 
valid for substructures. For example, the atoms bonded to the oxygen atom of ethanol are:

```@example structure
heavy = select(ethanol, "not element H")
O = heavy[findfirst(at -> name(at) == "O1", heavy)]
[ name(ethanol[j]) for j in ethanol.connectivity[index(O)] ]
```

When writing, `write_pdb` and `write_mmcif` number the atoms sequentially in the file (1, 2, ...)
and convert the connectivity to these numbers. Bonds to atoms that are not written are dropped:

```@example structure
tmpfile = tempname() * ".pdb"
write_pdb(tmpfile, heavy)
read_pdb(tmpfile).connectivity
```

!!! note
    The connectivity cannot be written if two of the atoms written have the same index (for 
    example, in the concatenation of two copies of a structure), because the bonds would be 
    ambiguous. In this case, renumber the atoms (and their connectivity), or skip writing the 
    connectivity with `connectivity=false`.

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
