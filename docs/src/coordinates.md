```@meta
CollapsedDocStrings = true
```

# Coordinate manipulation

## Atom positions

The `position`, `positions` and `set_positions!` functions are used to retrieve of set the coordinates of an atom (following the `AtomBase.jl` convention):

```@docs
position
set_position!
```

## Positions of arrays of Atoms

Use the `positions` function:

```@docs
positions(::AbstractVector{<:Atom})
```

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> positions(atoms[1:2])
2-element Vector{StaticArraysCore.SVector{3, Float32}}:
 [-9.229, -14.861, -5.481]
 [-10.048, -15.427, -5.569]
```

The `positions` function accepts selections:

C``\alpha`` coordinates:

```julia-repl
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> positions(atoms, "name CA")
3-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [-8.483, -14.912, -6.726]
 [-5.113, -13.737, -5.466]
 [-3.903, -11.262, -8.062]
```

The coordinates are output as arrays of static arrays (more specifically, as a `Vector{SVector{3,Float64}}`, from `StaticArrays`). 

## Move atoms and center of mass

The `center_of_mass` function can be used to compute the center of mass of set of atoms, and the 
`moveto!` function can be used to move the center of mass of the atoms to the origin (by default) 
or to a specified position:

```@docs
center_of_mass
moveto!
```

## Distance between sets of atoms

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

```@docs
distance
```

A special related function to compute the distance between a pair of residues is:

```@docs
residue_residue_distance
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

```@docs
closest
```

## Maximum and minimum coordinates of the atoms

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

```@docs
maxmin
```

## [Read and convert unit cells](@id read-unitcell)

```@docs
read_unitcell
lattice_to_matrix
matrix_to_lattice
```



