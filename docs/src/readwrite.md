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
PDBTools.Atom(1, 1, "N", "MET", "A", 1, 38.95, 49.3, 34.38, 0.0, 0.0, 1)

```

The data in the `Atom` structure is organized as follows:
```julia
struct Atom
  index :: Int64 # the sequential index of the atom in the file
  index_pdb :: Int64 # the number written in the index field of the PDB 
  name :: String
  resname :: String
  chain :: String
  resnum :: Int64
  x :: Float64
  y :: Float64
  z :: Float64
  b :: Float64
  occup :: Float64
  model :: Int64
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
 PDBTools.MutableAtom(1, 1, "N", "MET", "X", 1, 38.95, 49.3, 34.38, 0.0, 0.0, 1)
 PDBTools.MutableAtom(2, 2, "H1", "MET", "X", 1, 39.84, 49.01, 34.79, 0.0, 0.0, 1)
 PDBTools.MutableAtom(3, 3, "H2", "MET", "X", 1, 38.69, 48.52, 33.79, 0.0, 0.0, 1)
 PDBTools.MutableAtom(4, 4, "H3", "MET", "X", 1, 38.99, 50.19, 33.92, 0.0, 0.0, 1)
 ...

julia> atoms[1].resname = "AAA"
"AAA"

julia> atoms[1]
PDBTools.MutableAtom(1, 1, "N", "AAA", "X", 1, 38.95, 49.3, 34.38, 0.0, 0.0, 1)

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

