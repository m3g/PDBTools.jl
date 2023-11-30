# [Selection functions](@id selections)

A simple selection syntax is provided. Use it with, for example: 

```julia
atoms = select(atoms,"protein and resnum < 30")
```

## General selections 

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


### Retrieving indexes only 

If only the indexes of the atoms are of interest, a specific function
will directly return them:

```julia
indexes = selindex(atoms,"protein and name CA")

```

!!! note
    All indexing is 1-based. Thus, the first atom of the structure is atom 1.

## Use Julia functions directly

Selections can be done using Julia functions directly, providing a greater
control over the selection and, possibly, the use of user defined selection 
functions. For example:

```julia
myselection(atom) = (atom.x < 10.0 && atom.resname == "GLY") || (atom.name == "CA") 
atoms = select(atoms, by = myselection)
```
or, for example, using Julia anonymous functions
```julia
select(atoms, by = at -> isprotein(at) && name(at) == "O" && atom.x < 10.0)
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

The `eachresidue` iterator allows iteration over the resiudes of a structure (in PDB files distinct molecules are associated to different residues, thus this iterates similarly over the molecules of a structure). For example:

```julia
julia> protein = wget("1LBD")
   Array{Atoms,1} with 1870 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     SER     A      225        1   45.228   84.358   70.638 67.05  1.00     1       -         1
       2   CA     SER     A      225        1   46.080   83.165   70.327 68.73  1.00     1       -         2
       3    C     SER     A      225        1   45.257   81.872   70.236 67.90  1.00     1       -         3
                                                       ⋮ 
    1868  OG1     THR     A      462      238  -27.462   74.325   48.885 79.98  1.00     1       -      1868
    1869  CG2     THR     A      462      238  -27.063   71.965   49.222 78.62  1.00     1       -      1869
    1870  OXT     THR     A      462      238  -25.379   71.816   51.613 84.35  1.00     1       -      1870

julia> nALA = 0
       for residue in eachresidue(protein)
         if name(residue) == "ALA"
           nALA += 1
         end
       end
       nALA
22

```

The result of the iterator can also be collected, with:
```julia
julia> residues = collect(eachresidue(protein))
   Array{Residue,1} with 238 residues.


julia> residues[1]
 Residue of name SER with 6 atoms.
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     SER     A      225        1   45.228   84.358   70.638 67.05  1.00     1       -         1
       2   CA     SER     A      225        1   46.080   83.165   70.327 68.73  1.00     1       -         2
       3    C     SER     A      225        1   45.257   81.872   70.236 67.90  1.00     1       -         3
       4    O     SER     A      225        1   45.823   80.796   69.974 64.85  1.00     1       -         4
       5   CB     SER     A      225        1   47.147   82.980   71.413 70.79  1.00     1       -         5
       6   OG     SER     A      225        1   46.541   82.639   72.662 73.55  1.00     1       -         6

```

These residue vector *do not* copy the data from the original atom vector. Therefore, changes performed on these vectors will be reflected on the original data.  

It is possible also to iterate over the atoms of one or more residue:
```julia
julia> m_ALA = 0.
       for residue in eachresidue(protein)
         if name(residue) == "ALA"
           for atom in residue
             m_ALA += mass(atom)
           end
         end
       end
       m_ALA
1452.8601999999983

```
Which, in this simple example, results in the same as:

```julia
julia> sum( mass(atom) for atom in protein if atom.resname == "ALA" )
1452.8601999999983
```

## Using VMD

[VMD](https://www.ks.uiuc.edu/Research/vmd/) is a very popular and
powerful package for visualization of simulations. It contains a very
versatile library to read topologies and trajectory files, and a
powerful selection syntax. We provide here a wrapper to VMD which allows
using its capabilities.  

For example, the solute can be defined with: 
```julia
indexes, names = select_with_vmd("./system.gro","protein",vmd="/usr/bin/vmd")
solute = Selection(indexes,names,nmols=1)
```
The main advantage here is that all the file types that VMD supports are
supported. But VMD needs to be installed and is run in background, and
it takes a few seconds.     

The `select_with_vmd` function also accepts an optional keyword parameter `srcload`,
which can be used to load custom scripts within `vmd` before running setting
the selection. This allows the definition of `tcl` scripts with custom selection
macros, for instance. The usage would be: 
```julia
sel = select_with_vmd("file.pdb", "resname MYRES"; srcload = [ "mymacros1.tcl", "mymacros2.tcl" ])
```
Which corresponds to `source`ing each of the macro files in VMD before defining the 
selection with the custom `MYRES` name.

!!! compat
    The `select_with_vmd` function was introduced in PDBTools version 0.15.0.

!!! warning
    VMD uses 0-based indexing and `select_with_vmd` adjusts that. However, if
    a selection is performed by index, as with `index 1`, VMD will
    select the second atom, and the output will be `[2]`. Selections by
    type, name, segment, residue name, etc, will be consistent with one-based
    indexing.










