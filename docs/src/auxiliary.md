```@meta
CollapsedDocStrings = true
```

# Auxiliary functions 

## Get a protein sequence

To obtain a list of the residue names of the protein with three- and one-letter codes, use
```jldoctest
julia> using PDBTools, LazyArtifacts

julia> pdb = read_pdb(PDBTools.test_file("small.pdb"));

julia> getseq(pdb)
3-element Vector{String}:
 "A"
 "C"
 "D"
```

Use `getseq(atoms,code=2)` to get the sequence as three-letter residue codes, or `code=3` to get 
full natural-aminoacid names, like "Alanine", "Proline", etc:

```jldoctest
julia> using PDBTools, LazyArtifacts

julia> pdb = read_pdb(PDBTools.test_file("small.pdb"));

julia> getseq(pdb; code=2)
3-element Vector{String}:
 "ALA"
 "CYS"
 "ASP"

julia> getseq(pdb; code=3)
3-element Vector{String}:
 "Alanine"
 "Cysteine"
 "Aspartic acid"
```

```@docs
getseq
Sequence
```

!!! note
    If there is some non-standard protein residue in the sequence,
    inform the `getseq` function by adding a selection:
    ```jldoctest
    julia> using PDBTools

    julia> atoms = read_pdb(PDBTools.test_file("small.pdb"));

    julia> for at in atoms
              if resname(at) == "ALA"
                  at.resname = "NEW"
              end
           end

    julia> getseq(atoms, "protein or resname NEW"; code=2)
    3-element Vector{String}:
     "NEW"
     "CYS"
     "ASP"
    ```
    By default the selection will only return the sequence of natural amino acids. 

The `getseq` function can of course be used on an `Atom` list, accepts selections as the
last argument, as well as the reading and writing functions:

```jldoctest
julia> using PDBTools, LazyArtifacts

julia> atoms = read_pdb(PDBTools.test_file("small.pdb"));

julia> getseq(atoms, "residue > 1")
2-element Vector{String}:
 "C"
 "D"
```

## Residue tick labels for plots

The `residue_ticks` function provides a practical way to define tick labels in plots associated to an amino-acid sequence:

```julia
residue_ticks(
    atoms (or) residues (or) residue iterator; 
    first=nothing, last=nothing, stride=1, oneletter=true, serial=false,
)
```

The input structure can be provided as a vector of atoms (type `Vector{<:Atom}`) a residue iterator (obtained by `eachresidue(atoms)`) or a vector of residues (obtained by `collect(eachresidue(atoms))`). 

The function returns a tuple with residue numbers and residue names for the given atoms, to be used as tick labels in plots.

`first` and `last` optional keyword parameters are integers that refer to the residue numbers to be included. 
The `stride` option can be used to skip residues and declutter the tick labels.

If `oneletter` is `false`, three-letter residue codes are returned. Residues with unknown names will be 
named `X` or `XXX`. 

If `serial=false`, the positions of the ticks will be returned as the serial residue index in the structure.
If `serial=true`, the positions of the ticks are returned as their residue numbers. This difference is important
if the residue numbers do not start at `1`, and depending on the indexing of the data to be plotted.  

```@docs
residue_ticks
oneletter
threeletter
```

### Example

Here we illustrate how to plot the average temperature factor of each residue of a crystallographic model as function of the residues.

```@example residue_ticks
using PDBTools, LazyArtifacts, CairoMakie
atoms = read_mmcif(PDBTools.test_file("1UBQ.cif"), "protein")
residue_ticks(atoms; stride=10) # example of output
```

```@example residue_ticks
residues = collect(eachresidue(atoms))
fig = Figure(size = (600, 350))
ax = Axis(fig[1, 1];
    xlabel = "Residue",
    ylabel = "b-factor",
    xticks = residue_ticks(atoms; stride=10), # here we define the x-tick labels
    xticklabelrotation = π / 3,
)
lines!(ax, resnum.(residues), [sum(beta, r) / length(r) for r in residues]) # average b-factor per residue
fig
```

Alternatively (and sometimes conveniently), the residue ticks can be obtained by providing, 
instead of the `atoms` array, the residue iterator or the residue vector, as:

```@example residue_ticks
residue_ticks(eachresidue(atoms); stride=10)
```

```@example residue_ticks
residue_ticks(collect(eachresidue(atoms)); stride=10)
```

## Add hydrogens with OpenBabel

```@docs
add_hydrogens!
```

## Shared function names

PDBTools.jl is built on PDBToolsBase.jl, which defines the core data structures and the 
function names shared by the packages of the ecosystem (MolSimToolkit.jl, ComplexMixtures.jl, 
and others), which add methods to them.

```@docs
PDBTools.PDBToolsBase
PDBTools.distances
PDBTools.coordination_number
```

## Test files

The examples of this documentation use test files, which are downloaded on the first use:

```@docs
PDBTools.test_file
```
