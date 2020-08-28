# Atom properties

Some simple atom properties can be retrieved using special functions, which
operate on atoms of the type `Atom`. For example:

```julia
julia> atoms = readPDB("./file.pdb");

julia> atoms[1]
PDBTools.Atom(1, 1, "N", "ALA", "A", 1, -9.229, -14.861, -5.481, 0.0, 1.0, 1, 0)

julia> mass(atoms[1])
14.0067

julia> name(atoms[1])
"Nitrogen"

julia> atomic_number(atoms[1])
7

julia> element(atoms[1])
"N"

```
