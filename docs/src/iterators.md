```@meta
CollapsedDocStrings = true
```

# Iterators

PDBTools.jl provides lazy iterators over Residues, Segments, and Models of a structure file. The iterators behave similarly, and can be used bo computed properties of independent structural elements.

## Iterate over residues (or molecules)

The `eachresidue` iterator allows iteration over the residues of a structure (in PDB files distinct molecules are associated to different residues, thus this iterates similarly over the molecules of a structure). For example:

```jldoctest
julia> using PDBTools

julia> protein = read_pdb(PDBTools.SMALLPDB);

julia> count(atom -> resname(atom) == "ALA", protein)
12

julia> count(res -> resname(res) == "ALA", eachresidue(protein))
1
```

The result of the iterator can also be collected, with:
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

These residue vector *do not* copy the data from the original atom vector. Therefore, changes performed on these vectors will be reflected on the original data.  

It is possible also to iterate over the atoms of one or more residue:
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

Which, in this simple example, results in the same as:

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

## Iterate over segments 

The `eachsegment` iterator allows iteration over the segments of a structure. For example:

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
2-element Vector{Float64}:
 25222.339099999943
  1210.7300999999993

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

The `eachmodel` iterator allows iteration over the segments of a structure. For example:

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
11-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [0.633762128213737, -0.1413050285597195, -0.21796044955626692]
 [0.560772763043067, -0.15154922049365185, 0.1354801245061217]
 [0.506559232784597, -0.09771757024270422, 0.030405317843908077]
 ⋮
 [0.3889973654414868, -0.2110381926238272, 0.21802466991599198]
 [0.6995386823110438, -0.1537225338789714, 0.21793134264425737]

```

### Reference documentation

```@docs
Model
eachmodel
```

