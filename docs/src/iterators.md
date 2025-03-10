```@meta
CollapsedDocStrings = true
```

# Iterators

## Iterate over residues (or molecules)

```@docs
Residue
eachresidue
resname
residuename
```

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
```julia-repl
julia> using PDBTools

julia> protein = read_pdb(PDBTools.SMALLPDB);

julia> m_ALA = 0.
       for residue in eachresidue(protein)
         if name(residue) == "ALA"
           for atom in residue
             m_ALA += mass(atom)
           end
         end
       end
       m_ALA
73.09488999999999
```
Which, in this simple example, results in the same as:

```julia-repl
julia> sum(mass(at) for at in protein if resname(at) == "ALA" )
73.09488999999999
```

or

```julia-repl
julia> sum(mass(res) for res in eachresidue(protein) if resname(res) == "ALA" )
73.09488999999999
```

## Iterate over segments 

```@docs
Segment
eachsegment
```

The `eachsegment` iterator allows iteration over the segments of a structure. For example:

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> eachsegment(ats)
 Iterator with 2 segments.

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

julia> ats = read_pdb(PDBTools.DIMERPDB)

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
