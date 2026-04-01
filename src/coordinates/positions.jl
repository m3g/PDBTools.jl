import MolSimToolkitShared: positions

"""
    positions(atoms::AbstractVector{<:Atom}; selection) 

Returns the coordinates of the atoms. The input may be one atom (type `Atom`), a vector of atoms, or a `Residue`. 
The coordinates are returned as a vector of static vectors (from `StaticArrays`), more specifically
as a `Vector{SVector{3,Float64}}`.

### Examples

```julia-repl
julia> using PDBTools, StaticArrays 

julia> protein = wget("1LBD");

julia> positions(protein[1])
3-element SVector{3, Float64} with indices SOneTo(3):
 45.228
 84.358
 70.638

julia> positions(protein[1], as=SVector{3,Float32})
3-element SVector{3, Float32} with indices SOneTo(3):
 45.228
 84.358
 70.638

julia> positions(protein, "index <= 2")
2-element Vector{SVector{3, Float64}}:
 [45.228, 84.358, 70.638]
 [46.08, 83.165, 70.327]

julia> positions(protein, at -> at.resname == "ALA")
110-element Vector{SVector{3, Float64}}:
 [43.94, 81.982, 70.474]
 [43.02, 80.825, 70.455]
 [41.996, 80.878, 69.34]
 ⋮
 [-17.866, 84.088, 51.741]
 [-18.496, 83.942, 52.777]
 [-15.888, 82.583, 51.706]
  
julia> residues = collect(eachresidue(protein));

julia> positions(residues[1])
6-element Vector{SVector{3, Float64}}:
 [45.228, 84.358, 70.638]
 [46.08, 83.165, 70.327]
 [45.257, 81.872, 70.236]
 [45.823, 80.796, 69.974]
 [47.147, 82.98, 71.413]
 [46.541, 82.639, 72.662]

```

"""
positions(atoms::AbstractVector{<:Atom}) = position.(atoms)
function positions(atoms::AbstractVector{<:Atom}, selection_function::Function)
    x = SVector{3,Float32}[]
    for atom in atoms
        !selection_function(atom) && continue
        push!(x, position(atom))
    end
    return x
end
function positions(atoms::AbstractVector{<:Atom}, selection::AbstractString)
    return positions(atoms, parse_query(selection))
end


#
# Coordinates of the atoms of a residue/molecule
#
positions(residue::Residue, selection_function::Function=all) = positions(residue.atoms[residue.range], selection_function)
positions(residue::Residue, selection::AbstractString) = positions(residue.atoms[residue.range], selection)

@testitem "positions" begin
    atoms = read_pdb(PDBTools.TESTPDB)
    s = select(atoms, "residue = 3")
    @test all(isapprox.(stack(positions(s)), [
            -4.383 -4.51 -3.903 -3.731 -4.938 -4.417 -5.543 -5.867 -5.451 -6.974 -2.626 -1.94
            -11.903 -11.263 -11.262 -12.076 -10.279 -9.552 -9.911 -10.85 -10.837 -11.289 -10.48 -10.014
            -6.849 -6.096 -8.062 -8.767 -8.612 -9.06 -7.784 -9.684 -10.863 -9.3 -7.749 -8.658
        ]; atol=1e-3))
    r = Residue(select(atoms, "residue = 3"))
    @test positions(s) == positions(r)
    residues = collect(eachresidue(atoms))
    @test positions(select(atoms, "residue = 3")) == positions(residues[3])
    @test positions(atoms, "residue = 3") == positions(s)
    @test positions(residues[1]) == positions(select(atoms, "residue = 1"))
    @test all(isapprox.(positions(residues[1], at -> name(at) == "N")[1], [-9.229, -14.861, -5.481]; atol=1e-3))
    @test all(isapprox.(positions(residues[1], "name N")[1], [-9.229, -14.861, -5.481]; atol=1e-3))
end


