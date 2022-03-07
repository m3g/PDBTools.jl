"""

```
coor(atoms; selection) 
```

Returns the coordinates of the atoms. The input may be one atom (type `Atom`), a vector of atoms, or a `Residue`. 
The coordinates are returned as a vector of static vectors (from `StaticArrays`), more specifically
as a `Vector{SVector{3,Float64}}`.

### Examples

```julia-repl
julia> using PDBTools, StaticArrays 

julia> protein = wget("1LBD");

julia> coor(protein[1])
3-element SVector{3, Float64} with indices SOneTo(3):
 45.228
 84.358
 70.638

julia> coor(protein[1],as=SVector{3,Float32})
3-element SVector{3, Float32} with indices SOneTo(3):
 45.228
 84.358
 70.638

julia> coor(protein, "index <= 2")
2-element Vector{SVector{3, Float64}}:
 [45.228, 84.358, 70.638]
 [46.08, 83.165, 70.327]

julia> coor(protein, only = at -> at.resname == "ALA")
110-element Vector{SVector{3, Float64}}:
 [43.94, 81.982, 70.474]
 [43.02, 80.825, 70.455]
 [41.996, 80.878, 69.34]
 ⋮
 [-17.866, 84.088, 51.741]
 [-18.496, 83.942, 52.777]
 [-15.888, 82.583, 51.706]
  
julia> residues = collect(eachresidue(protein));

julia> coor(residues[1])
6-element Vector{SVector{3, Float64}}:
 [45.228, 84.358, 70.638]
 [46.08, 83.165, 70.327]
 [45.257, 81.872, 70.236]
 [45.823, 80.796, 69.974]
 [47.147, 82.98, 71.413]
 [46.541, 82.639, 72.662]

```

"""
coor(atom::Atom) = SVector{3,Float64}(atom.x, atom.y, atom.z)

function coor(atoms::AbstractVector{Atom}, selection::String)
    query = parse_query(selection)
    return coor(atoms, only = atom -> apply_query(query, atom))
end

function coor(atoms::AbstractVector{Atom}; only = all)
    n = 0
    for atom in atoms
        !only(atom) && continue
        n = n + 1
    end
    x = fill(zero(SVector{3,Float64}),n)
    i = 0
    for atom in atoms
        !only(atom) && continue
        i += 1
        x[i] = coor(atom)
    end
    return x
end

#
# Coordinates of the atoms of a residue/molecule
#
coor(residue::Residue; only = all) = coor(residue.atoms[residue.range], only = only)
coor(residue::Residue, selection::String) = coor(residue.atoms[residue.range], selection)
