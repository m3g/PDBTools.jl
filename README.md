# PDBTools
Simple structure and functions to read and write PDB files

## Installing:

```
julia> ] add https://github.com/mcubeg/PDBTools
```

## Fundamental structure: Atom

`Atom` is the structure of data containing the atom index, name,
residue, coordinates, etc.

## Functions provided:

`readPDB(filename)`: Reads a PDB file and returns a vector of atoms of
type `Atom`.

`editPDB(filename)`: Reads a PDB file and returns a vector of atoms
of type `MutableAtom`, which contains the same data, but the data can be
modified.

`read_atom(pdb_line)`: Given a line of a PDB file containing atom data,
returns the data in a `Atom` structure.

`write_atom(atom::Atom)`: Given an atom in the `Atom` structure, returns
a string formatted in the PDB format, to be written to a file. 




