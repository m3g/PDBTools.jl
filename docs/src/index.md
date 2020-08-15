# PDBTools

PDBTools is a simple package to read and write Protein Data Bank files,
select atoms, and work with their coordinates.  

### Features:

> Simple data structure: 
> ```julia
>  Atom(5, 5, "CA", "ALA", "A", 1, -8.483, -14.912, -6.726, 0.0, 1.0, 1, 0)
> ```

> Selection syntax:
> ```julia
> resname ARG and name CA
> ```

> Allows use of Julia (possibly user-defined) functions for selection:
> ```julia
> atom -> ( atom.resname == "ARG" && atom.x < 10 ) || atom.name == "N"
> ```

### Not indicated for:

We do not aim to provide the fastest PDB parsing methods. If
speed in reading files, returning subsets of the structures, etc., is
critical to you, probably you will do better with some packages of 
[BioJulia](https://github.com/BioJulia), 
[BioStructures](https://github.com/BioJulia/BioStructures.jl) in
particular.

