# PDBTools
Simple structures and functions to read and write PDB files

## Installing:

```
julia> ] add https://github.com/m3g/PDBTools
```

## Using:

```
julia> using PDBTools
```

## Fundamental structure: Atom

`Atom` is the structure of data containing the atom index, name,
residue, coordinates, etc.

## Functions provided:

`PDBTools.readPDB(filename)`: Reads a PDB file and returns a vector of atoms of
type `Atom`.

`PDBTools.editPDB(filename)`: Reads a PDB file and returns a vector of atoms
of type `MutableAtom`, which contains the same data, but the data can be
modified.

`PDBTools.writePDB(atoms,filename)`: Writes a PDB file of name `filename` from the vector
`atoms` which contains the atom data in the `Atom` or `MutableAtom` structures.

`PDBTools.read_atom(pdb_line)`: Given a line of a PDB file containing atom data,
returns the data in a `Atom` structure.

`PDBTools.write_atom(atom::Atom)`: Given an atom in the `Atom` structure, returns
a string formatted in the PDB format, to be written to a file. 

`PDBTools.getseq(filename or Vector{Atom})`: Returns a list of residue names with three
and one letter codes.

# Selection functions

`PDBTools.select(atoms,"protein and resnum < 30")`: Simple selection tool that provides the
most common selection features. Here, `atoms` is a vector obtained with `readPDB` or `editPDB`.

`PDBTools.selindex(atoms,"protein and resnum < 30")`: Simple selection tool that provides the
most common selection features. Here, `atoms` is a vector obtained with `readPDB` or `editPDB`.
This one returns only the atomic indexes of corresponding to the selection.
