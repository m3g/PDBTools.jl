# [Selection functions](@id selections)

A simple selection syntax is provided. Use it with, for example: 

```julia
atoms = PDBTools.select(atoms,"protein and resnum < 30")
```

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
| `protein`  |                       |             | `protein`     |
| `water`    |                       |             | `water`       |
| `beta`     | `=`,`>`,`<`,`<=`,`>=` | Real        | `beta > 0.5`  |
| `occup`    | `=`,`>`,`<`,`<=`,`>=` | Real        | `occup >= 0.3`|


### Retrieving indexes only 

If only the indexes of the atoms are of interest, a specific function
will directly return them:

```julia
indexes = PDBTools.selindex(atoms,"protein and name CA")

```

!!! note
    All indexing is 1-based. Thus, the first atom of the structure is atom 1.

