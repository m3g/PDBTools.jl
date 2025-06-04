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

## Using VMD

[VMD](https://www.ks.uiuc.edu/Research/vmd/) is a very popular and
powerful package for visualization of simulations. It contains a very
versatile library to read topologies and trajectory files, and a
powerful selection syntax. We provide here a wrapper to VMD which enables
using its capabilities.  

The `select_with_vmd` input can be a vector of `PDBTools.Atom`s, or a filename.
If the input is a vector of `Atom`s, the output will be the corresponding atoms
matching the selection. If the input is a filename, two lists are returned:
the list of indices and names of the corresponding atoms. This is because some
input files supported by VMD (e. g. `GRO`, `PSF`, etc.) do not contain full atom
information. 

For example, here some atoms are selected from a previously loaded vector of atoms:

```julia-repl
julia> using PDBTools

julia> pdbfile = PDBTools.SMALLPDB

julia> atoms = read_pdb(pdbfile);

julia> selected_atoms = select_with_vmd(atoms,"resname ALA and name HT2 HT3";vmd="/usr/bin/vmd")
   Vector{Atom{Nothing}} with 2 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
       4  HT3     ALA     A        1        1   -8.652  -15.208   -4.741  0.00  0.00     1    PROT         4
```

And, now, we provide the filename as input:

```julia-repl
julia> selected_atoms = select_with_vmd(pdbfile,"resname ALA and name HT2 HT3";vmd="/usr/bin/vmd")
([3, 4], ["HT2", "HT3"])
```

Note that in the above examples we use `name HT2 HT3` which is not currently supported by the 
internal PDBTools `select` function, which would require `name HT2 or name HT3`. 

Here, the output will contain two lists, one of atom indices (one-based) and atom names.
The indices correspond to sequential indices in the input, *not* the indices written in the PDB file, for example.

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
the selection. This enables the definition of `tcl` scripts with custom selection
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










