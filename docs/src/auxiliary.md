```@meta
CollapsedDocStrings = true
```

# Some auxiliary functions to quickly retrieve some data 

## Get the protein sequence

```@docs
getseq
```

To obtain a list of the residue names of the protein with three- and one-letter codes, use
```jldoctest
julia> using PDBTools

julia> getseq(PDBTools.SMALLPDB)
3-element Vector{String}:
 "A"
 "C"
 "D"
```

Use `getseq(atoms,code=2)` to get the sequence as three-letter residue codes, or `code=3` to get 
full natural-aminoacid names, like "Alanine", "Proline", etc:

```jldoctest
julia> using PDBTools

julia> getseq(PDBTools.SMALLPDB, code=2)
3-element Vector{String}:
 "ALA"
 "CYS"
 "ASP"

julia> getseq(PDBTools.SMALLPDB, code=3)
3-element Vector{String}:
 "Alanine"
 "Cysteine"
 "Aspartic acid"
```

!!! note
    If there is some non-standard protein residue in the sequence,
    inform the `getseq` function by adding a selection:
    ```jldoctest
    julia> using PDBTools

    julia> atoms = read_pdb(PDBTools.SMALLPDB);

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
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> getseq(atoms, "residue > 1")
2-element Vector{String}:
 "C"
 "D"
```

## Distance between sets of atoms

```@docs
distance
```

The distance between atoms, or sets of atoms, can be computed with the `distance` function. This
function returns the *minimum distance* between the atoms of the sets involved. For example:

```julia-repl
julia> using PDBTools

julia> model = wget("1BSX");

julia> protein = select(model,"protein");

julia> ligand = select(model,"resname T3");

julia> distance(protein,ligand)
2.7775834820937417
```

## Closest atoms and their distance

A function similar to the one above is `closest`, which returns the shortest distance between atoms
but also the identity of the atom or pair of atoms that satisfy that shortest distance:

```julia-repl
julia> using PDBTools

julia> model = wget("1BSX");

julia> protein = select(model,"protein");

julia> ligand = select(model,"resname T3");

julia> closest(ligand,protein)
(43, 3684, 2.7775834820937417)

julia> ligand[43]
    4037   O1      T3     B        2      512  -22.568   81.625    3.159  1.00 36.59     1       -      4041

julia> protein[3684]
    3684  NE2     HIS     B      435      472  -21.539   82.145    5.686  1.00 44.44     1       -      3686

julia> distance(ligand[43],protein[3684])
2.7775834820937417
```

## Obtain arrays with coordinates

```@docs
coor
```

Use the `coor` function:

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> coor(atoms[1])
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  -9.229
 -14.861
  -5.481

julia> coor(atoms[1:2])
2-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [-9.229, -14.861, -5.481]
 [-10.048, -15.427, -5.569]
```

The `coor` function accepts selections:

C``\alpha`` coordinates:

```julia-repl
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> coor(atoms, "name CA")
3-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [-8.483, -14.912, -6.726]
 [-5.113, -13.737, -5.466]
 [-3.903, -11.262, -8.062]
```

The coordinates are output as arrays of static arrays (more specifically, as a `Vector{SVector{3,Float64}}`, from `StaticArrays`). 

## Maximum and minimum coordinates of the atoms

```@docs
maxmin
```

Use `maxmin(atoms)`, or `maxmin(atoms,"resname CA")`, for example:

```julia-repl
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> maxmin(atoms, "residue > 1")
 Minimum atom coordinates: xmin = [-6.974, -16.785, -10.863]
 Maximum atom coordinates: xmax = [-1.94, -9.552, -3.844]
 Length in each direction: xlength = [5.034000000000001, 7.2330000000000005, 7.019]
```

`m` is a structure containing the three vectors with minimum and maximum
coordinates, and lengths.

## Residue tick labels for plots

```@docs
residue_ticks
```

The `residue_ticks` function provides a practical way to define tick labels in plots associated to an amino-acid sequence:

```julia
residue_ticks(
    atoms (or) residues (or) residue iterator; 
    first=nothing, last=nothing, stride=1, oneletter=true, serial=false,
)
```

The input structure can be provided as a vector of atoms (type `Vector{Atom}`) a residue iterator (obtained by `eachresidue(atoms)`) or a vector of residues (obtained by `collect(eachresidue(atoms))`). 

The function returns a tuple with residue numbers and residue names for the given atoms, to be used as tick labels in plots.

`first` and `last` optional keyword parameters are integers that refer to the residue numbers to be included. 
The `stride` option can be used to skip residues and declutter the tick labels.

If `oneletter` is `false`, three-letter residue codes are returned. Residues with unknown names will be 
named `X` or `XXX`. 

if `serial=false` the positions of the ticks will be returned as a the serial residue index in the structure.
If `serial=true` the positions of the ticks are returned as their residue numbers. This difference is important
if the residue numbers do not start at `1` and depending on the indexing of the data to be plotted.  

!!! compat
    The functionality of the `residue_ticks` as described requires PDBTools version 1.6.0 or greater. 

    The `serial` option was introduced in v1.8.1

### Example

Here we illustrate how to plot the average temperature factor of each residue of a crystallographic model as function of the residues.

```julia-repl
julia> using PDBTools, Plots

julia> atoms = wget("1UBQ", "protein");

julia> residue_ticks(atoms; stride=10) # example of output
([1, 11, 21, 31, 41, 51, 61, 71], ["M1", "K11", "D21", "Q31", "Q41", "E51", "I61", "L71"])

julia> plot(
           resnum.(eachresidue(atoms)), # x-axis: residue numbers
           [ mean(beta.(res)) for res in eachresidue(atoms) ], # y-axis: average b-factor per residue
           xlabel="Residue", 
           xticks=residue_ticks(atoms; stride=10), # here we define the x-tick labels
           ylabel="b-factor", 
           xrotation=60,
           label=nothing, framestyle=:box,
      )
```

Produces the following plot:

![./assets/residue_ticks.png](./assets/residue_ticks.png)

Alternatively (and sometimes conveniently), the residue ticks can be obtained by providing, 
instead of the `atoms` array, the residue iterator or the residue vector, as:

```julia-repl
julia> residue_ticks(eachresidue(atoms); stride=10)
([1, 11, 21, 31, 41, 51, 61, 71], ["M1", "K11", "D21", "Q31", "Q41", "E51", "I61", "L71"])

julia> residue_ticks(collect(eachresidue(atoms)); stride=10)
([1, 11, 21, 31, 41, 51, 61, 71], ["M1", "K11", "D21", "Q31", "Q41", "E51", "I61", "L71"])
```
