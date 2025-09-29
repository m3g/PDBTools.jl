```@meta
CollapsedDocStrings = true
```

# Iterators

PDBTools.jl provides lazy iterators over Residues, Chains, Segments, and Models of a structure file. The iterators behave similarly, and can be used bo computed properties of independent structural elements.
The documentation bellow exemplifies in more detail the features associated to Residue and Chain interators, but the properties and valid for Segment and Model iterators similarly.

## Iterate over residues (or molecules)

The `eachresidue` iterator enables iteration over the residues of a structure. In PDB files, distinct molecules are often treated as separate residues, so this iterator can be used to iterate over the molecules within a structure. For example:

```jldoctest
julia> using PDBTools

julia> protein = read_pdb(PDBTools.SMALLPDB);

julia> count(atom -> resname(atom) == "ALA", protein)
12

julia> count(res -> resname(res) == "ALA", eachresidue(protein))
1
```

Here, the first `count` counts the number of atoms with the residue name "ALA", while the second uses `eachresidue` to count the number of residues named "ALA". This highlights the distinction between residue-level and atom-level operations.

### Collecting Residues into a Vector

Residues produced by `eachresidue` can be collected into a vector for further processing:

```jldoctest
julia> using PDBTools

julia> protein = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(protein))
3-element Vector{Residue}[
    ALA1A
    CYS2A
    ASP3A
]

julia> residues[1]
 Residue of name ALA with 12 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
       2 1HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
⋮
      11    C     ALA     A        1        1   -7.227  -14.047   -6.599  1.00  0.00     1    PROT        11
      12    O     ALA     A        1        1   -7.083  -13.048   -7.303  1.00  0.00     1    PROT        12
```

!!! note
    Iterators or collected vectors *do not* create copies of the original atom data. This means that any changes made to the residue vector will directly modify the corresponding data in the original atom vector.

### Iterating Over Atoms Within Residues

You can iterate over the atoms of one or more residues using nested loops. Here, we compute the total number of atoms of ALA residues: 

```jldoctest
julia> using PDBTools

julia> protein = read_pdb(PDBTools.SMALLPDB);

julia> n_ala_cys = 0
       for residue in eachresidue(protein)
            if name(residue) in ("ALA", "CYS")
                for atom in residue
                   n_ala_cys += 1
                end
            end
       end
       n_ala_cys
23
```

This method produces the same result as the more concise approach:

```jldoctest 
julia> using PDBTools

julia> protein = read_pdb(PDBTools.SMALLPDB);

julia> sum(length(r) for r in eachresidue(protein) if name(r) in ("ALA", "CYS"))
23
```

### Reference documentation

```@docs
Residue
eachresidue
resname
residuename
```

## Iterate over chains

The `eachchain` iterator in PDBTools enables users to iterate over the chains in a PDB structure. A PDB file may contain multiple protein chains. This iterator simplifies operations involving individual chains.


```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.CHAINSPDB);

julia> chain.(eachchain(ats)) # Retrieve the names of all chains in the structure
4-element Vector{InlineStrings.String3}:
 "A"
 "B"
 "A"
 "D"

julia> model.(eachchain(ats)) # Retrieve the model numbers associated with each chain
4-element Vector{Int32}:
 1
 1
 1
 2

julia> chain_A1 = first(eachchain(ats)); # Access the first chain in the iterator

julia> resname.(eachresidue(chain_A1)) # Retrieve residue names for chain A in model 1
3-element Vector{InlineStrings.String7}:
 "ASP"
 "GLN"
 "LEU"

```
In the example above, the `chain.` command retrieves the names of all chains in the structure, while  `model.` command lists the model numbers for each chain. This PDB structure contains two models for chain A, where the third residue changes from leucine (LEU) in model 1 to valine (VAL) in model 2.

### Accessing Chains by Index

As seen in the previous example, The `first` and `last` commands allow quick access to the first an last elements in the iterator. For more specific indexing, you can collect all chains into an array and then use numerical indices to access them.

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.CHAINSPDB);

julia> chains = collect(eachchain(ats))
4-element Vector{Chain}[
    Chain(A-48 atoms)
    Chain(B-48 atoms)
    Chain(A-48 atoms)
    Chain(D-45 atoms)
]

julia> chain_B = chains[2]
 Chain B with 48 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
      49    N     ASP     B        4        4  135.661  123.866  -22.311  1.00  0.00     1    ASYN        49
      50   CA     ASP     B        4        4  136.539  123.410  -21.227  1.00  0.00     1    ASYN        50
⋮
      95 HD23     LEU     B        6        6  138.780  120.216  -17.864  1.00  0.00     1    ASYN        95
      96    O     LEU     B        6        6  141.411  117.975  -21.923  1.00  0.00     1    ASYN        96

```

### Modifying Atom Properties in a Chain

Any changes made to the atoms of a chain variable directly overwrite the properties of the original atoms in the structure. For example, modifying the occupancy and beta-factor columns of atoms in model 2 of chain A will update the corresponding properties in the original structure.

In the example below, the `occup` and `beta` properties of all atoms in model 2 of chain A are set to 0.00. The changes are reflected in the original `ats` vector, demonstrating that the modifications propagate to the parent data structure.

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.CHAINSPDB);

julia> first(eachchain(ats))
 Chain A with 48 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ASP     A        1        1  133.978  119.386  -23.646  1.00  0.00     1    ASYN         1
       2   CA     ASP     A        1        1  134.755  118.916  -22.497  1.00  0.00     1    ASYN         2
⋮
      47 HD23     LEU     A        3        3  130.568  111.868  -26.242  1.00  0.00     1    ASYN        47
      48    O     LEU     A        3        3  132.066  112.711  -21.739  1.00  0.00     1    ASYN        48

 
julia> for chain in eachchain(ats)
           if name(chain) == "A" && model(chain) == 2
               for atom in chain
                   atom.occup = 0.00
                   atom.beta = 0.00
               end
           end
       end

julia> first(eachchain(ats))
 Chain A with 48 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ASP     A        1        1  133.978  119.386  -23.646  1.00  0.00     1    ASYN         1
       2   CA     ASP     A        1        1  134.755  118.916  -22.497  1.00  0.00     1    ASYN         2
⋮
      47 HD23     LEU     A        3        3  130.568  111.868  -26.242  1.00  0.00     1    ASYN        47
      48    O     LEU     A        3        3  132.066  112.711  -21.739  1.00  0.00     1    ASYN        48

```

This behavior ensures efficient data manipulation but requires careful handling to avoid unintended changes. 

### Reference documentation

```@docs
Chain
eachchain
```

## Iterate over segments 

The `eachsegment` iterator enables iteration over the segments of a structure. For example:

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> eachsegment(ats)
 Segment iterator with length = 2

julia> name.(eachsegment(ats))
2-element Vector{InlineStrings.String7}:
 "A"
 "B"
```

The result of the iterator can also be collected, with:
```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> s = collect(eachsegment(ats))
2-element Vector{Segment}[ 
    A-(1905 atoms))
    B-(92 atoms))
]

julia> s[1]
 Segment of name A with 1905 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     LYS     A      211        1   52.884   24.022   35.587  1.00 53.10     1       A         1
       2   CA     LYS     A      211        1   52.916   24.598   36.993  1.00 53.10     1       A         2
⋮
    1904  OD2     ASP     A      461      243   17.538   51.009   45.748  1.00 97.43     1       A      1904
    1905  OXT     ASP     A      461      243   14.506   47.082   47.528  1.00 97.43     1       A      1905
```

These segment structure *does not* copy the data from the original atom vector. Therefore, changes performed on these vectors will be reflected on the original data.  

Iterators can be used to obtain or modify properties of the segments. Here we illustrate computing the mass of
each segment and renaming segment of all atoms with the segment indices:

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> s = collect(eachsegment(ats))
2-element Vector{Segment}[ 
    A-(1905 atoms))
    B-(92 atoms))
]

julia> mass.(s)
2-element Vector{Float32}:
 25222.553
  1210.7296

julia> for (iseg, seg) in enumerate(eachsegment(ats))
           for at in seg
               at.segname = "$(at.segname)$iseg"
           end
       end

julia> collect(eachsegment(ats))
2-element Vector{Segment}[ 
    A1-(1905 atoms))
    B2-(92 atoms))
]
```

### Reference documentation

```@docs
Segment
eachsegment
```

## Iterate over models

The `eachmodel` iterator enables iteration over the segments of a structure. For example:

```jldoctest
julia> using PDBTools

julia> ats = wget("8S8N");

julia> eachmodel(ats)
 Model iterator with length = 11

julia> model.(eachmodel(ats))
11-element Vector{Int32}:
  1
  2
  3
  ⋮
 10
 11
```

The result of the iterator can also be collected, with:

```jldoctest
julia> using PDBTools

julia> ats = wget("8S8N");

julia> m = collect(eachmodel(ats))
11-element Vector{Model}[
    1-(234 atoms))
    2-(234 atoms))
    ⋮
    10-(234 atoms))
    11-(234 atoms))
]

julia> m[1]
 Model 1 with 234 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     DLE     A        2        1   -5.811   -0.380   -2.159  1.00  0.00     1                 1
       2   CA     DLE     A        2        1   -4.785   -0.493   -3.227  1.00  0.00     1                 2
⋮
     233  HT2   A1H5T     B      101       13   -5.695    5.959   -3.901  1.00  0.00     1               233
     234  HT1   A1H5T     B      101       13   -4.693    4.974   -2.743  1.00  0.00     1               234
```

The model structure *does not* copy the data from the original atom vector. Therefore, changes performed on these vectors will be reflected on the original data.  

Iterators can be used to obtain or modify properties of the segments. Here we illustrate computing the mass of
each segment and renaming segment of all atoms with the segment indices:

```jldoctest
julia> using PDBTools

julia> ats = wget("8S8N");

julia> center_of_mass.(eachmodel(ats))
11-element Vector{StaticArraysCore.SVector{3, Float32}}:
 [0.6337627, -0.14130484, -0.2179606]
 [0.56077266, -0.15154965, 0.1354806]
 [0.5065595, -0.0977174, 0.030405657]
 ⋮
 [0.38899764, -0.21103837, 0.2180245]
 [0.69953984, -0.15372278, 0.21793146]

```

### Reference documentation

```@docs
Model
eachmodel
```

