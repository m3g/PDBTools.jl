# [Selection functions](@id selections)

A simple selection syntax is provided. Use it with, for example: 

```julia
atoms = PDBTools.select(atoms,"protein and resnum < 30")
```

## General selections 

Accepted boolean operators: `and`, `or`, and `not`. 

The accepted keywords for the selection are: 

| Keyword    | Options               | Input value | Example       | 
|:----------:|:---------------------:|:-----------:|:-------------:|
| `index`    | `=`,`>`,`<`,`<=`,`>=` | Integer     | `index <= 10` |
| `name`     |                       | String      | `name CA`     |
| `resname`  |                       | String      | `resname ALA` |
| `resnum`   | `=`,`>`,`<`,`<=`,`>=` | Integer     | `resnum = 10` |
| `chain`    |                       | String      | `chain A`     |
| `model`    |                       | Integer     | `model 1`     |
| `b`        | `=`,`>`,`<`,`<=`,`>=` | Real        | `b > 0.5`     |
| `occup`    | `=`,`>`,`<`,`<=`,`>=` | Real        | `occup >= 0.3`|

## Special macros: proteins, water

Just use these keywords to select the residues mathching the properties
desired. 

Examples:
```julia
aromatic = PDBTools.select("file.pdb","aromatic")

```
```julia
atoms = PDBTools.select("file.pdb") 
aromatic = PDBTools.select(atoms,"charged")

```

Available keywords:

| Keyword       |
|:-------------:|
| `water`       |
| `protein`     |
| `backbone`    |
| `sidechain`   |
| `acidic`      |
| `aliphatic`   |
| `aromatic`    |
| `basic`       |
| `charged`     |
| `hydrophobic` |
| `neutral`     |
| `nonpolar`    |
| `polar`       |

!!! note  
    The properties refer to protein residues and will return `false`
    to every non-protein residue. Thus, be careful with the use of `not`
    with these selections, as they might retrieve non-protein atoms.


### Retrieving indexes only 

If only the indexes of the atoms are of interest, a specific function
will directly return them:

```julia
indexes = PDBTools.selindex(atoms,"protein and name CA")

```

!!! note
    All indexing is 1-based. Thus, the first atom of the structure is atom 1.

