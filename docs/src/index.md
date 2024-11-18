# PDBTools

PDBTools is a simple package to read and write Protein Data Bank files,
select atoms, and work with their coordinates. It is aimed to provide support
for the typical uses of structure files in the context of molecular dynamics
simulations.

As of version 2.0, PDBTools is able to read and write the atomic data 
from PDB and mmCIF structure files.

## Features:

Simple data structure: 
```julia-repl
julia> printatom(atoms[1])
  index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
      1   OW     SOL     X        1        1   54.370   45.310   33.970  0.00  0.00     1       -         1
```

Selection syntax:
```julia
resname ARG and name CA
```

Allows use of Julia (possibly user-defined) functions for selection:
```julia
atom -> ( atom.resname == "ARG" && atom.x < 10 ) || atom.name == "N"
```
### Not indicated for:

PDBTools is not very strict in following the PDB or mmCIF formats. In particular,
it does not read any of the meta-data of these files, only `ATOM` and `HETATM` fields
are of interest. Also, it supports repeated atom entries, as each atom is read as 
an independent object. This flexibility provides support for common structure formats
occurring in the Molecular Dynamics Simulations field. 

If more comprehensive (and strict) support for these files is necessary, use the packages of 
[BioJulia](https://github.com/BioJulia), 
[BioStructures](https://github.com/BioJulia/BioStructures.jl) in
particular.

