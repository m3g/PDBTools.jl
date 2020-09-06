# PDBTools
Simple structures and functions to read and write PDB files

## Documentation

The documentation can be found at: [http://m3g.iqm.unicamp.br/PDBTools](http://m3g.iqm.unicamp.br/PDBTools)

## Installing:

```
julia> ] add PDBTools

julia> using PDBTools

```
### Features:

> Simple data structure: 
> ```julia
>atoms = readPDB("./structure.pdb")
>   Array{PDBTools.Atom,1} with 62026 atoms with fields:
>   index name resname chain   resnum  residue        x        y        z     b occup model segname index_pdb
>       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
>       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
>       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
>                                                       ⋮ 
>   62023   H2    TIP3     C     9338    19637  -18.014   16.666   11.615  0.00  1.00     1    WAT2     62023
>   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  0.00  1.00     1    WAT2     62024
>   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025
>   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026
> 
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



