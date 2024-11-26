```@meta
CollapsedDocStrings = true
```

# [Selection functions](@id selections)

The `select` function can be used to select subsets of atoms from a vector
of atoms. A simple selection syntax can be used, for example: 

```julia
atoms = select(atoms, "protein and resnum < 30")
```

or standard Julia function can be provided as the second argument:

```julia
atoms = select(atoms, at -> isprotein(at) && resnum(at) < 30)
```

## General selection syntax 

Accepted Boolean operators: `and`, `or`, and `not`. 

The accepted keywords for the selection are: 

| Keyword    | Options               | Input value | Example       | 
|:----------:|:---------------------:|:-----------:|:-------------:|
| `index`    | `=`,`>`,`<`,`<=`,`>=` | Integer     | `index <= 10` |
| `index_pdb`| `=`,`>`,`<`,`<=`,`>=` | Integer     | `index_pdb <= 10` |
| `name`     |                       | String      | `name CA`     |
| `element`  |                       | String      | `element N`   |
| `resname`  |                       | String      | `resname ALA` |
| `resnum`   | `=`,`>`,`<`,`<=`,`>=` | Integer     | `resnum = 10` |
| `residue`  | `=`,`>`,`<`,`<=`,`>=` | Integer     | `residue = 10`|
| `chain`    |                       | String      | `chain A`     |
| `model`    |                       | Integer     | `model 1`     |
| `beta`     | `=`,`>`,`<`,`<=`,`>=` | Real        | `beta > 0.5`  |
| `occup`    | `=`,`>`,`<`,`<=`,`>=` | Real        | `occup >= 0.3`|
| `segname`  |                       | String      | `segname PROT`|
|            |                       |             |               |

!!! note
    `resnum` is the residue number as written in the PDB file, while `residue`
    is the residue number counted sequentially in the file.

    `index_pdb` is the number written in the "atom index" field of the PDB file,
    while `index` is the sequential index of the atom in the file. 


## Special macros: proteins, water

Just use these keywords to select the residues matching the properties
desired. 

Examples:
```julia
aromatic = select(atoms,"aromatic")

```
```julia
aromatic = select(atoms,"charged")

```

Available keywords:

| Keywords      |               |               |
|:-------------:|:-------------:|:-------------:|
| `water`       |               |               |
| `protein`     | `backbone`    | `sidechain`   |
| `acidic`      | `basic`       |               |
| `aliphatic`   | `aromatic`    |               |
| `charged`     | `neutral`     |               |
| `polar`       | `nonpolar`    |               |
| `hydrophobic` |               |               |
|               |               |               |

!!! note  
    The properties refer to protein residues and will return `false`
    to every non-protein residue. Thus, be careful with the use of `not`
    with these selections, as they might retrieve non-protein atoms.

```@docs
select
Select
```

## Retrieving indices, filtering, etc

If only the indices of the atoms are of interest, the Julia `findall`
function can be used, by passing a `Select` object, or a regular 
function, to select the atoms:

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.TESTPDB, "protein and residue <= 3");

julia> findall(Select("name CA"), atoms)
3-element Vector{Int64}:
  5
 15
 26

julia> findall(at -> name(at) == "CA", atoms)
3-element Vector{Int64}:
  5
 15
 26
```

!!! note
    All indexing is 1-based. Thus, the first atom of the structure is atom 1.

The `Select` constructor can be used to feed simple selection syntax entries to 
other Julia functions, such as `findfirst`, `findlast`, or `filter`:

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.TESTPDB, "protein and residue <= 3");

julia> filter(Select("name CA"), atoms)
   Vector{Atom{Nothing}} with 3 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       5   CA     ALA     A        1        1   -8.483  -14.912   -6.726  1.00  0.00     1    PROT         5
      15   CA     CYS     A        2        2   -5.113  -13.737   -5.466  1.00  0.00     1    PROT        15
      26   CA     ASP     A        3        3   -3.903  -11.262   -8.062  1.00  0.00     1    PROT        26

julia> findfirst(Select("beta = 0.00"), atoms)
1
```

!!! tip
    The `sel""` literal string macro is a shortcut for `Select`. Thus, these syntaxes are valid:
    ```jldoctest
    julia> using PDBTools

    julia> atoms = read_pdb(PDBTools.TESTPDB, "protein and residue <= 3");

    julia> name.(filter(sel"name CA", atoms))
    3-element Vector{InlineStrings.String7}:
     "CA"
     "CA"
     "CA"

    julia> findfirst(sel"name CA", atoms)
    5
    ```

## Use Julia functions directly

Selections can be done using Julia functions directly, providing a greater
control over the selection and, possibly, the use of user defined selection 
functions. For example:

```julia
myselection(atom) = (atom.x < 10.0 && atom.resname == "GLY") || (atom.name == "CA") 
atoms = select(atoms, myselection)
```
or, for example, using Julia anonymous functions
```julia
select(atoms, at -> isprotein(at) && name(at) == "O" && atom.x < 10.0)
```

The only requirement is that the function defining the selection receives an `PDBTools.Atom` as
input, and returns `true` or `false` depending on the conditions required for the atom.

!!! note
    The macro-keywords described in the previous section can be used within 
    the Julia function syntax, but the function names start with `is`. For example:
    ```julia
    select(atoms, at -> isprotein(at) && resnum(at) in [ 1, 5, 7 ])
    ```
    Thus, the macro selection functions are:
    `iswater`, 
    `isprotein`,     `isbackbone`,    `issidechain`,
    `isacidic`,      `isbasic`,                  
    `isaliphatic`,   `isaromatic`,               
    `ischarged`,     `isneutral`,                
    `ispolar`,       `isnonpolar`,               
    and `ishydrophobic`.                          
    


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
   Vector{Residue} with 3 residues.

julia> residues[1]
 Residue of name ALA with 12 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
       2 1HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
                                                       ⋮
      10  HB3     ALA     A        1        1   -9.164  -15.063   -8.765  1.00  0.00     1    PROT        10
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

```@docs
Residue
eachresidue
resname
residuename
```

## Iterate over chains

The `eachchain` iterator allows for iteration over the chains of a structure. A PDB file may contain more than one protein monomer listed as different chains, thus this iterator behaves similarly to iterating over each monomer of a structure. For example:


```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.CHAINSPDB);

julia> chain.(eachchain(ats))
3-element Vector{InlineStrings.String3}:
 "A"
 "B"
 "C"

```
The command `chain` provides the name of all the chains stored by the iterator `eachchain`. In this example, the PDB loaded has 3 chains identified in alphabetical order with the letters A, B and C. 

The chains can be stored in diferent variables. In the example below, all the atoms belonging to chain A are stored in the variable chain_A with the command `first`, which selects the first element of the iterator `eachchain`. The command `last` and the indexing list is also supported. However, this type of indexing can only be used in a vector of elements `Atom`, so the chains must first to be collected in an array. 

```julia-repl
julia> using PDBTools

ats = read_pdb(PDBTools.CHAINSPDB);

julia> chain_A = first(eachchain(ats))
 Chain of name A with 48 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ASP     A        1        1  133.978  119.386  -23.646  1.00  0.00     1    ASYN         1
       2   CA     ASP     A        1        1  134.755  118.916  -22.497  1.00  0.00     1    ASYN         2
       3    C     ASP     A        1        1  135.099  117.439  -22.652  1.00  0.00     1    ASYN         3
                                                       ⋮ 
      46 HD22     LEU     A        3        3  130.704  113.003  -27.586  1.00  0.00     1    ASYN        46
      47 HD23     LEU     A        3        3  130.568  111.868  -26.242  1.00  0.00     1    ASYN        47
      48    O     LEU     A        3        3  132.066  112.711  -21.739  1.00  0.00     1    ASYN        48

julia> chains = collect(eachchain(ats))
   Array{Chain,1} with 3 chains.

julia> chain_B = chains[2]
 Chain of name B with 48 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
      49    N     ASP     B        4        4  135.661  123.866  -22.311  1.00  0.00     1    ASYN        49
      50   CA     ASP     B        4        4  136.539  123.410  -21.227  1.00  0.00     1    ASYN        50
      51    C     ASP     B        4        4  137.875  122.934  -21.788  1.00  0.00     1    ASYN        51
                                                       ⋮ 
      94 HD22     LEU     B        6        6  139.485  119.501  -16.418  1.00  0.00     1    ASYN        94
      95 HD23     LEU     B        6        6  138.780  120.216  -17.864  1.00  0.00     1    ASYN        95
      96    O     LEU     B        6        6  141.411  117.975  -21.923  1.00  0.00     1    ASYN        96

julia> chain_C = last(chains)
 Chain of name C with 48 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
      97    N     ASP     C        7        7  137.110  128.213  -20.946  1.00  0.00     1    ASYN        97
      98   CA     ASP     C        7        7  137.938  127.668  -19.879  1.00  0.00     1    ASYN        98
      99    C     ASP     C        7        7  137.485  128.196  -18.521  1.00  0.00     1    ASYN        99
                                                       ⋮ 
     142 HD22     LEU     C        9        9  141.370  132.570  -11.550  1.00  0.00     1    ASYN       142
     143 HD23     LEU     C        9        9  142.275  131.839  -10.225  1.00  0.00     1    ASYN       143
     144    O     LEU     C        9        9  141.311  127.197  -12.279  1.00  0.00     1    ASYN       144


```

Any changes made to the atoms of a chain variable will overwrite the properties of the original atoms. In the example below, the occupancy column of the atoms in chain_A is set to the value of 0.00. When calling the original ats vector, it is possible to see that the occupancy column values of the same set of atoms changed from 1.00 to 0.00.

```julia-repl
julia> using PDBTools

julia> ats = read_pdb(PDBTools.CHAINSPDB)
   Vector{Atom{Nothing}} with 144 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ASP     A        1        1  133.978  119.386  -23.646  1.00  0.00     1    ASYN         1
       2   CA     ASP     A        1        1  134.755  118.916  -22.497  1.00  0.00     1    ASYN         2
       3    C     ASP     A        1        1  135.099  117.439  -22.652  1.00  0.00     1    ASYN         3
                                                       ⋮ 
     142 HD22     LEU     C        9        9  141.370  132.570  -11.550  1.00  0.00     1    ASYN       142
     143 HD23     LEU     C        9        9  142.275  131.839  -10.225  1.00  0.00     1    ASYN       143
     144    O     LEU     C        9        9  141.311  127.197  -12.279  1.00  0.00     1    ASYN       144

julia> chain_A = first(eachchain(ats));

julia> for atom in chain_A
       atom.occup = 0.00
       end

julia> ats
   Vector{Atom{Nothing}} with 144 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ASP     A        1        1  133.978  119.386  -23.646  0.00  0.00     1    ASYN         1
       2   CA     ASP     A        1        1  134.755  118.916  -22.497  0.00  0.00     1    ASYN         2
       3    C     ASP     A        1        1  135.099  117.439  -22.652  0.00  0.00     1    ASYN         3
                                                       ⋮ 
     142 HD22     LEU     C        9        9  141.370  132.570  -11.550  1.00  0.00     1    ASYN       142
     143 HD23     LEU     C        9        9  142.275  131.839  -10.225  1.00  0.00     1    ASYN       143
     144    O     LEU     C        9        9  141.311  127.197  -12.279  1.00  0.00     1    ASYN       144 

```

When iterating over chains, it is possible to make equal changes along every chain. For instance, each chain in the example below starts with the  aspartic acid residue, ASP, but if one wishes to equally rename only the first ASP in each chain to its protonated state name, ASPP - as it is implemented in the CHARMM22 force field -, one would do as follows:


```julia-repl
julia> using PDBTools

ats = read_pdb(PDBTools.CHAINSPDB);

julia> for chx in eachchain(ats)
       for res in first(eachresidue(chx))
       res.resname = "ASPP"
       end
       end

julia> count(res -> resname(res) == "ASPP", eachresidue(ats))
3 

```

```@docs
Chain
eachchain
```

## Using VMD

[VMD](https://www.ks.uiuc.edu/Research/vmd/) is a very popular and
powerful package for visualization of simulations. It contains a very
versatile library to read topologies and trajectory files, and a
powerful selection syntax. We provide here a wrapper to VMD which allows
using its capabilities.  

For example, the solute can be defined with: 
```julia
indices, names = select_with_vmd("./system.gro","protein",vmd="/usr/bin/vmd")
```
The output will contain two lists, one of atom indices (one-based) and atom names.
The indices correspond to sequential indices in the input, *not* the indices
written in the PDB file, for example.

The input may also be a vector of atoms of type `PDBTools.Atom`:

```julia
atoms = read_pdb("mypdbfile.pdb")
indices, names = select_with_vmd(atoms,"protein",vmd="/usr/bin/vmd")
```

!!! tip
    If `vmd` is available in your path, there is no need to pass it as a keyword parameter.

The main advantage here is that all the file types and the complete selection syntax 
that VMD supports are supported. But VMD needs to be installed and is run in background, and
it takes a few seconds to run.

```@docs
select_with_vmd
```

###  Loading vmd scripts

The `select_with_vmd` function also accepts an optional keyword parameter `srcload`,
which can be used to load custom scripts within `vmd` before running setting
the selection. This allows the definition of `tcl` scripts with custom selection
macros, for instance. The usage would be: 
```julia
sel = select_with_vmd("file.pdb", "resname MYRES"; srcload = [ "mymacros1.tcl", "mymacros2.tcl" ])
```
Which corresponds to `source`ing each of the macro files in VMD before defining the 
selection with the custom `MYRES` name.

!!! warning
    VMD uses 0-based indexing and `select_with_vmd` adjusts that. However, if
    a selection is performed by index, as with `index 1`, VMD will
    select the second atom, and the output will be `[2]`. Selections by
    type, name, segment, residue name, etc, will be consistent with one-based
    indexing.










