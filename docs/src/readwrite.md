# Read and write files

## Read a PDB file

To read a PDB file and return a vector of atoms of
type `Atom`, do:
```julia
atoms = readPDB("file.pdb")

```

`Atom` is the structure of data containing the atom index, name,
residue, coordinates, etc. For example, after reading a file (as shown
bellow), a list of atoms with the following structure will be generated:

```julia
julia> atoms[1]
PDBTools.Atom(1, 1, "N", "ALA", "A", 1, 1, -9.229, -14.861, -5.481, 0.0, 1.0, 1, "PROT")

```

The data in the `Atom` structure is organized as follows:
```julia
struct Atom
  index :: Int64 # The sequential index of the atoms in the file
  index_pdb :: Int64 # The index as written in the PDB file (might be anything)
  name :: String
  resname :: String
  chain :: String
  resnum :: Int64 # Number of residue as written in PDB file
  residue :: Int64 # Sequential residue (molecule) number in file
  x :: Float64
  y :: Float64
  z :: Float64
  b :: Float64
  occup :: Float64
  model :: Int64
  segname :: String # Segment name (cols 73:76)
end
```

!!! tip
    For all these reading and writting function, a final argument can be provided
    to read or write a subset of the atoms, following the selection syntax described 
    in the [Selection](@ref selections) section. For example:
    ```julia
    protein = readPDB("file.pdb","protein")

    ```
    or
    ```julia
    arginines = readPDB("file.pdb","resname ARG")

    ```
    The only difference is that, if using Julia anonymous functions, the
    keyword is `only`:
    ```julia
    arginines = readPDB("file.pdb",only = atom -> atom.resname == "ARG")

    ```
    The same is valid for `edit` and `write` functions, below. 
      

## Edit a PDB file

Using the `editPDB` function, a vector of the same structure as above is
returned, but of `MutableAtom` type, meaning that the content of every
field can be modified. For example:
```julia
julia> atoms = editPDB("file.pdb")
1500-element Array{PDBTools.Atom,1}:
 PDBTools.MutableAtom(1, 1, "N", "ALA", "A", 1, 1, -9.229, -14.861, -5.481, 0.0, 1.0, 1, "PROT")
 PDBTools.MutableAtom(2, 2, "HT1", "ALA", "A", 1, 1, -10.048, -15.427, -5.569, 0.0, 0.0, 1, "PROT")
 PDBTools.MutableAtom(3, 3, "HT2", "ALA", "A", 1, 1, -9.488, -13.913, -5.295, 0.0, 0.0, 1, "PROT")
 PDBTools.MutableAtom(4, 4, "HT3", "ALA", "A", 1, 1, -8.652, -15.208, -4.741, 0.0, 0.0, 1, "PROT")
 ...

julia> atoms[1].segname = "ABCD"
"AAA"

julia> atoms[1]
 PDBTools.MutableAtom(1, 1, "N", "ALA", "A", 1, 1, -9.229, -14.861, -5.481, 0.0, 1.0, 1, "ABCD")

```

## Write a PDB file

To write a PDB file use the `writePDB` function, as:

```julia
writePDB(atoms,"file.pdb")

```
where `atoms` contain a list of atoms in the `Atom` or `MutableAtom` structures.

# Read and write single-atom lines 

`PDBTools.read_atom(pdb_line)`: Given a line of a PDB file containing atom data,
returns the data in a `Atom` structure. To convert the `atom` read with
this function into a mutable structure, use `atom = MutableAtom(atom)`.

`PDBTools.write_atom(atom::Atom)`: Given an atom in the `Atom` structure, returns
a string formatted in the PDB format, to be written to a file. 

