# Some auxiliary functions to quickly retrive some data 

## Get the protein sequence

To obtain a list of the residue names of the protein with three- and one-letter codes, use
```julia
julia> seq = PDBTools.getseq("file.pdb")
76×2 Array{String,2}:
 "VAL"  "V"
 "LYS"  "K"
 ⋮      
 "ARG"  "R"
 "GLY"  "G"

```

!!! note
    If there is some non-standard protein residue in the sequence,
    inform the `getseq` function by adding a selection:
    ```julia
    julia> PDBTools.getseq("file.pdb","protein or resname NEW")
    76×2 Array{String,2}:
     "VAL"  "V"
     "NEW"  "N"
     ⋮      
     "ARG"  "R"
     "GLY"  "G"
    ```
    A list of new residues can also be provided, with `newres=["NEW1","NEW2",...]`

The `getseq` function can of course be used on an `Atom` list, accepts selections as the
last argument, as well as the reading and writting functions:

```julia
atoms = PDBTools.readPDB("file.pdb")
seq = PDBTools.getseq(atoms,"chain A")

```

## Obtain arrays with coordinates

All atoms:

```julia
julia> x = PDBTools.coor(atoms)
1476×3 Array{Float64,2}:
 38.03  49.56  35.45
 38.12  52.85  37.52
  ⋮
 60.05  47.5   57.34
 63.52  46.9   58.93

```

Or use selections to retrive the coordinates of subsets of atoms:

C``\alpha`` coordinates:

```julia
julia> xCA = PDBTools.coor(atoms,"name CA")
76×3 Array{Float64,2}:
 38.03  49.56  35.45
 38.12  52.85  37.52
  ⋮
 60.05  47.5   57.34
 63.52  46.9   58.93

```
