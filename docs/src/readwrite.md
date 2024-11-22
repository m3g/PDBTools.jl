```@meta
CollapsedDocStrings = true
```

# Read and write files

PDBTools can read and write `PDB` and `mmCIF` files. The relevant functions are:

```@docs
read_pdb
read_mmcif
```

!!! note 
    In the following examples, the `read_pdb` function will be illustrated. The usage is
    similar to that of `read_mmcif`, to read `mmCIF (PDBx)` files. 

## Read a PDB/mmCIF file

To read a PDB file and return a vector of atoms of
type `Atom`, do:
```julia
atoms = read_pdb("file.pdb")
```

`Atom` is the structure of data containing the atom index, name,
residue, coordinates, etc. For example, after reading a file (as shown
bellow), a list of atoms with the following structure will be generated:

```julia-repl
julia> printatom(atoms[1])
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
```

The data in the `Atom` structure is organized as indicated in the following documentation:

```@docs
Atom
```

!!! tip
    For all these reading and writing functions, a final argument can be provided
    to read or write a subset of the atoms, following the selection syntax described 
    in the [Selection](@ref selections) section. For example:
    ```julia
    protein = read_pdb("file.pdb","protein")
    ```
    or
    ```julia
    arginines = read_pdb("file.pdb","resname ARG")
    ```
    The only difference is that, if using Julia anonymous functions, the
    keyword is `only`:
    ```julia
    arginines = read_pdb("file.pdb", only = atom -> atom.resname == "ARG")
    ```
    The same is valid for the `write` function, below. 

## [Atom field assignment in mmCIF files](@id field_assignment)

By default, the assignment of the `_atom_site` fields of the mmCIF format to the fields of the `Atom` data structure 
follows the [standard mmCIF convention](https://mmcif.wwpdb.org/docs/tutorials/content/atomic-description.html):

```julia
Dict{String,Symbol}(
    "id" => :index_pdb
    "Cartn_x" => :x
    "Cartn_y" => :y
    "Cartn_z" => :z
    "occupancy" => :occup
    "B_iso_or_equiv" => :beta
    "pdbx_formal_charge" => :charge
    "pdbx_PDB_model_num" => :model
    "label_atom_id" => :name
    "label_comp_id" => :resname
    "label_asym_id" => :chain
    "auth_seq_id" => :resnum
    "type_symbol" => :pdb_element
)
```

This assignment can be customized by providing the `field_assignment` keyword parameter to the `read_mmcif` function. 
In the following example, we exemplify the possibility of reading `_atom_site.type_symbol` field of the mmCIF file into the `name` field of the
atom data structure:

```jldoctest
julia> using PDBTools

julia> ats = read_mmcif(PDBTools.TESTCIF, "index <= 5");

julia> name.(ats)
5-element Vector{InlineStrings.String7}:
 "N"
 "CA"
 "C"
 "O"
 "N"

julia> ats = read_mmcif(PDBTools.TESTCIF, "index <= 5"; 
           field_assignment=Dict("type_symbol" => :name)
        );

julia> name.(ats)
5-element Vector{InlineStrings.String7}:
 "N"
 "C"
 "C"
 "O"
 "N"

```

The custom entries set in the `field_assignment` keyword will overwrite the default 
assignments for entries sharing keys or fields. For instance, in the example above,
the `label_atom_id` fields which is by default assigned to `:name` is not being read
anymore.

## Get structure from the Protein Data Bank

Use the `wget` function to retrieve the atom data directly from the PDB database,
optionally filtering the atoms with a selection:

```julia-repl
julia> atoms = wget("1LBD","name CA")
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       2   CA     SER     A      225        1   46.080   83.165   70.327 68.73  1.00     1       -         2
       8   CA     ALA     A      226        2   43.020   80.825   70.455 63.69  1.00     1       -         8
      13   CA     ASN     A      227        3   41.052   82.178   67.504 53.45  1.00     1       -        13
                                                       ⋮
    1847   CA     GLN     A      460      236  -22.650   79.082   50.023 71.46  1.00     1       -      1847
    1856   CA     MET     A      461      237  -25.561   77.191   51.710 78.41  1.00     1       -      1856
    1864   CA     THR     A      462      238  -26.915   73.645   51.198 82.96  1.00     1       -      1864
```

```@docs
wget
```

## Write a PDB/mmCIF file

To write a PDB file use the `write_pdb` function, as:

```julia
write_pdb("file.pdb", atoms)
```
where `atoms` contain a list of atoms with the `Atom` structures.

```@docs
write_pdb
write_mmcif
```

The use of the `field_assignment` keyword, as explained in the [field assignment](@ref field_assignment) section
is possible in the call to `write_mmcif`. 

# Read from string buffer

In some cases a PDB file data may be available as a string and not a regular file. For example,
when reading the output of a zipped file. In these cases, it is possible to obtain the array
of atoms by reading directly the string buffer with, for example:

```jldoctest
julia> using PDBTools

julia> pdbdata = read(PDBTools.TESTPDB, String); # returns a string with the PDB data, to exemplify

julia> atoms = read_pdb(IOBuffer(pdbdata), "protein and name CA")
   Vector{Atom{Nothing}} with 104 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       5   CA     ALA     A        1        1   -8.483  -14.912   -6.726  1.00  0.00     1    PROT         5
      15   CA     CYS     A        2        2   -5.113  -13.737   -5.466  1.00  0.00     1    PROT        15
      26   CA     ASP     A        3        3   -3.903  -11.262   -8.062  1.00  0.00     1    PROT        26
                                                       ⋮
    1425   CA     GLU     A      102      102    4.414   -4.302   -7.734  1.00  0.00     1    PROT      1425
    1440   CA     CYS     A      103      103    4.134   -7.811   -6.344  1.00  0.00     1    PROT      1440
    1454   CA     THR     A      104      104    3.244  -10.715   -8.603  1.00  0.00     1    PROT      1454
```

## Edit a Vector{<:Atom} object

The `Atom` structure is mutable, meaning that the fields can be edited. For example:

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.TESTPDB)
   Vector{Atom{Nothing}} with 62026 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
                                                       ⋮
   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  1.00  0.00     1    WAT2     62024
   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  1.00  0.00     1    WAT2     62025
   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  1.00  0.00     1    WAT2     62026

julia> atoms[1].segname = "ABCD"
"ABCD"

julia> printatom(atoms[1])
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    ABCD         1
```

Additionally, With the `edit!` function, you can directly edit or view the data in a
vector of `Atoms` in your preferred text editor. 

```julia-repl
julia> edit!(atoms)
```

This will open a text editor. Here, we modified the data in the `resname` field of the first atom
to `ABC`. Saving and closing the file will update the `atoms` array:

```julia-repl
julia> printatom(atoms[1])
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ABC     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
```

```@docs
edit!
```



