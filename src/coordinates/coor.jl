"""
    coor(atoms; selection) 

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
coor(atom::Atom) = SVector{3,Float32}(atom.x, atom.y, atom.z)

function coor(atoms::AbstractVector{<:Atom}, selection::String)
    query = parse_query(selection)
    return coor(atoms, only=atom -> apply_query(query, atom))
end

function coor(atoms::AbstractVector{<:Atom}; only=all)
    x = SVector{3,Float32}[]
    for atom in atoms
        !only(atom) && continue
        push!(x, coor(atom))
    end
    return x
end

#
# Coordinates of the atoms of a residue/molecule
#
coor(residue::Residue; only=all) = coor(residue.atoms[residue.range], only=only)
coor(residue::Residue, selection::String) = coor(residue.atoms[residue.range], selection)

@testitem "coor" begin
    atoms = read_pdb(PDBTools.TESTPDB)
    s = select(atoms, "residue = 3")
    @test all(isapprox.(stack(coor(s)), [
            -4.383 -4.51 -3.903 -3.731 -4.938 -4.417 -5.543 -5.867 -5.451 -6.974 -2.626 -1.94
            -11.903 -11.263 -11.262 -12.076 -10.279 -9.552 -9.911 -10.85 -10.837 -11.289 -10.48 -10.014
            -6.849 -6.096 -8.062 -8.767 -8.612 -9.06 -7.784 -9.684 -10.863 -9.3 -7.749 -8.658
        ]; atol=1e-3))
    r = Residue(select(atoms, "residue = 3"))
    @test coor(s) == coor(r)
    residues = collect(eachresidue(atoms))
    @test coor(select(atoms, "residue = 3")) == coor(residues[3])
    @test coor(atoms, "residue = 3") == coor(s)
    @test coor(residues[1]) == coor(select(atoms, "residue = 1"))
    @test all(isapprox.(coor(residues[1]; only=at -> name(at) == "N")[1], [-9.229, -14.861, -5.481]; atol=1e-3))
    @test all(isapprox.(coor(residues[1], "name N")[1], [-9.229, -14.861, -5.481]; atol=1e-3))
end
