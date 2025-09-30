import CellListMap
using CellListMap: ParticleSystem, map_pairwise!
using LinearAlgebra
using Statistics: mean

# Container for the custom atom filed that will carry the atom SASA
struct SASA
    sasa::Float32
end

#=
    generate_dots(atomi_radius, probe_radius, n_dots)

Generates points on the surface of a sphere of a given radius (atomic_radius + probe_radius) 
using the double cubic lattice method. `n_dots` controls the maximum number of dots
per sphere.

=#
function generate_dots(atomic_radius, probe_radius::Real, n_dots::Int)
    radius = atomic_radius + probe_radius
    n_grid_points = round(Int, (3 / (4 * π * ((1 / 2)^3)) * n_dots)^(1 / 3))
    if radius <= 0 || n_grid_points <= 0
        throw(ArgumentError("""\n
            Too small or incorrect probe_radius ($probe_radius) or number of dots ($n_dots).

        """))
    end

    dots = SVector{3,Float32}[]
    step = 2 * radius / n_grid_points
    radius_sq = radius^2

    # --- First Lattice ---
    for i in 0:n_grid_points, j in 0:n_grid_points, k in 0:n_grid_points
        x = -radius + i * step
        y = -radius + j * step
        z = -radius + k * step
        dist_sq = x^2 + y^2 + z^2

        # Project internal points onto the surface
        if 0 < dist_sq < radius_sq
            norm_factor = radius / sqrt(dist_sq)
            push!(dots, norm_factor * SVector{3,Float32}(x, y, z))
        end
    end

    # --- Second Lattice (shifted by half a step) ---
    offset = step / 2
    for i in 0:n_grid_points, j in 0:n_grid_points, k in 0:n_grid_points
        x = -radius + offset + i * step
        y = -radius + offset + j * step
        z = -radius + offset + k * step
        dist_sq = x^2 + y^2 + z^2

        # Project internal points onto the surface
        if dist_sq < radius_sq
            norm_factor = radius / sqrt(dist_sq)
            push!(dots, norm_factor * SVector{3,Float32}(x, y, z))
        end
    end

    # Return unique points to avoid redundancy
    return unique!(dots)
end

#
# Structure that carries the information about dots of each atom, if 
# they are exposed or found to be occluded by other atoms.
#
struct AtomDots
    exposed::BitVector
end
function CellListMap.reset_output!(x::AtomDots)
    x.exposed .= true
    return x
end
CellListMap.copy_output(x::AtomDots) = AtomDots(copy(x.exposed))
function CellListMap.reducer(x::AtomDots, y::AtomDots)
    x.exposed .= x.exposed .& y.exposed
    return x
end

function update_dot_exposure!(x, y, dot_cache_i, surface_dots_i, rj_sq)
    for idot in eachindex(surface_dots_i.exposed)
        if surface_dots_i.exposed[idot]
            dot_on_surface = x + dot_cache_i[idot]
            # Position the dot on the atom's surface in the molecule's coordinate system
            # Check if the dot is inside the neighboring atom j
            if sum(abs2, dot_on_surface - y) < rj_sq
                surface_dots_i.exposed[idot] = false
            end
        end
    end
    return nothing
end

function update_pair_dot_exposure!(
    x, y, i, j, d2, surface_dots;
    atoms, dot_cache, atom_type::F1, atom_radius_from_type::F2, probe_radius,
) where {F1,F2}
    type_i = atom_type(atoms[i])
    type_j = atom_type(atoms[j])
    ri = atom_radius_from_type(type_i) + probe_radius
    rj = atom_radius_from_type(type_j) + probe_radius
    ri_sq = ri^2
    rj_sq = rj^2
    if d2 < (ri_sq + rj_sq + 2*ri*rj)
        dot_cache_i = dot_cache[type_i]
        dot_cache_j = dot_cache[type_j]
        surface_dots_i = surface_dots[i]
        surface_dots_j = surface_dots[j]
        update_dot_exposure!(x, y, dot_cache_i, surface_dots_i, rj_sq)
        update_dot_exposure!(y, x, dot_cache_j, surface_dots_j, ri_sq)
    end
    return surface_dots
end

"""
    atomic_sasa(atoms; probe_radius, n_dots)

Calculates the Solvent Accessible Surface Area (SASA) for a vector of `Atom`s. 

# Main argument

- `atoms::Vector{PDBTools.Atom}`: A vector of atoms in the molecule.

# Returns

Vector of `PDBToools.Atom`s, with custom fields of type `PDBTools.SASA`, 
which contain the solvent accessible surface area for each atom, in Å².
The `sasa` function computes the total SASA or the SASA of a subset of the atoms
in the structure.

# Optional arguments 

- `probe_radius::Real=1.4`: The radius of the solvent probe in Angstroms.
- `n_dots::Int=100`: The number of grid points along one axis for dot generation. 
  Higher values lead to more accurate but slower calculations.
- `parallel::Bool=true`: Control if the computation runs in parallel (requires 
  running Julia with multiple threads).

# Example

```jldoctest
julia> using PDBTools

julia> prot = select(read_pdb(PDBTools.TESTPDB), "protein");

julia> at_sasa = atomic_sasa(prot);

julia> sasa(at_sasa) # total sasa of prot
5352.8564f0

julia> sasa(at_sasa, "backbone") # backbone sasa in prot
982.7195f0

julia> sasa(at_sasa, "not backbone") # other atoms
4370.136f0

julia> sasa(at_sasa, "resname ARG GLU") # some residue types
546.56476f0
```

# Additional control:

Two arguments can control the atom radii used for computing the SASA. These arguments
are functions:

- `atom_type`: Function that given each atom of the array of atoms, returns the atom "type".
- `atom_radius_from_type`: Given the atom "type", returns the vdW radius of the atom. 

By default, `atom_type = PDBTools.element`, a function that just returns the element symbol of the atom,
and `atom_radius_from_type` obtains the vdW radius from the `PDBTools.elements` list given the element symbol.
Here, the atomc radii of https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page) are used. 
Atoms with missing radius have a `NaN` value, and the computation will not return meaningful
values. 

*Algorithm:* Eisenhaber F, Lijnzaad P, Argos P, Sander C, & Scharf M (1995) J. Comput. Chem. 16, 273-284.

"""
function atomic_sasa(
    atoms::AbstractVector{<:Atom};
    probe_radius::Real=1.4,
    n_dots::Int=100,
    atom_type::Function=element,
    atom_radius_from_type::Function=type -> getproperty(elements[type], :vdw_radius),
    parallel=true,
)

    # Unique list of atom types
    atom_types = atom_type.(unique(atom_type, atoms))

    # Memoization for dot generation to avoid recomputing for same radii
    dot_cache = Dict{eltype(atom_types),Vector{SVector{3,Float32}}}()
    for type in atom_types
        atom_radius = atom_radius_from_type(type)
        if isnan(atom_radius)
            throw(ArgumentError("""\n
                Atom of type $type does not have a vdW radius defined.
                Use custom `atom_type` and `atom_radius_from_type` input parameters if needed. 
            """))
        end
        dot_cache[type] = generate_dots(atom_radius_from_type(type), probe_radius, n_dots)
    end

    system = ParticleSystem(
        xpositions=coor.(atoms),
        unitcell=nothing,
        cutoff=2 * (maximum(atom_radius_from_type(type) for type in atom_types) + probe_radius),
        output=[AtomDots(trues(length(dot_cache[atom_type(at)]))) for at in atoms],
        output_name=:surface_dots,
        parallel=parallel,
    )

    map_pairwise!(
        (x, y, i, j, d2, surface_dots) ->
            update_pair_dot_exposure!(
                x, y, i, j, d2, surface_dots;
                atoms, dot_cache, atom_type, atom_radius_from_type, probe_radius,
            ),
        system,
    )

    ats_with_sasa = [
        add_custom_field(
            atoms[i],
            SASA(
                4 * π * (atom_radius_from_type(atom_type(atoms[i])) + probe_radius)^2 *
                mean(system.surface_dots[i].exposed)
            )
        ) for i in eachindex(atoms)
    ]

    return ats_with_sasa
end

"""
    sasa(atom::Atom{PDBTools.SASA})
    sasa(atoms::AbstractVector{PDBTools.SASA})
    sasa(atoms::AbstractVector{PDBTools.SASA}, selection::Union{Function,String})

Given the output of `atomic_sasa`, sums up contributions of atoms to compute the SASA
of the full structure, an atom, or a subset of atoms. The function can called with a 
single `Atom{SASA}` atom, a vector of atoms (in which case the full SASA is returned),
or a vector of atoms and a selection, given by a function or selection string. 

*Algorithm:* Eisenhaber F, Lijnzaad P, Argos P, Sander C, & Scharf M (1995) J. Comput. Chem. 16, 273-284

# Example

```jldoctest
julia> using PDBTools

julia> prot = select(read_pdb(PDBTools.TESTPDB), "protein");

julia> at_sasa = atomic_sasa(prot);

julia> sasa(at_sasa) # total sasa of prot
5352.8564f0

julia> sasa(at_sasa, "backbone") # selection string
982.7195f0

julia> sasa(at_sasa, at -> name(at) == "CA") # selection function
41.121826f0

julia> sasa(at_sasa[1]) # single atom SASA
5.806664f0
```

"""
function sasa end

sasa(atom::Atom{SASA}) = atom.custom.sasa
sasa(atoms::AbstractVector{<:Atom{SASA}}) = sum(sasa, atoms)
sasa(atoms::AbstractVector{<:Atom{SASA}}, sel::String) = sasa(atoms, Select(sel))
sasa(atoms::AbstractVector{<:Atom{SASA}}, sel::Function) = sum(sasa(at) for at in atoms if sel(at))

@testitem "sasa" begin

    using PDBTools
    prot = read_pdb(PDBTools.TESTPDB, "protein")

    # Compare with the output of VMD: the difference is that VMD uses a vdW radius for H of 1.00, and we use 1.10.
    const vmd_radii = Dict(
        "N" => 1.55,
        "C" => 1.70,
        "H" => 1.00,
        "O" => 1.52,
        "S" => 1.80,
    )
    at_sasa = atomic_sasa(prot; atom_radius_from_type=type -> vmd_radii[type])
    @test sasa(at_sasa) ≈ 5364.78515625 rtol = 0.1
    # Accessiblity of groups within the structure
    @test sasa(at_sasa, "backbone") ≈ 1128.9261474609375 rtol = 0.05
    @test sasa(at_sasa, "resname GLU LYS") ≈ 796.6279296875 rtol = 0.05
    @test sasa(at_sasa, "residue 1") ≈ 124.51239776611328 rtol = 0.05
    @test sasa(at_sasa, "residue 104") ≈ 122.6167221069336 rtol = 0.05

    # Compare with Gromacs output
    at_sasa = atomic_sasa(prot)
    # Isolated groups
    @test sasa(atomic_sasa(select(prot, "backbone and not name O"))) ≈ 100 * 55.245 rtol = 0.05
    @test sasa(atomic_sasa(select(prot, "name CA"))) ≈ 100 * 58.642 rtol = 0.05
    @test sasa(atomic_sasa(select(prot, "sidechain and not element H"))) ≈ 100 * 69.029 rtol = 0.05

    # Test non-contiguous indexing with general selections
    at_sasa = atomic_sasa(select(prot, "name CA"))
    @test sasa(at_sasa) ≈ 5856.476f0
    @test sasa(at_sasa, "resname THR") ≈ 325.33087f0
    @test sasa(at_sasa, at -> resname(at) == "THR") ≈ 325.33087f0

    # Test parallelization
    @test sasa(atomic_sasa(prot; parallel=false)) ≈ sasa(atomic_sasa(prot; parallel=true))

    # Test errors
    @test_throws ArgumentError atomic_sasa(prot; probe_radius=-2.0)
    @test_throws ArgumentError atomic_sasa(prot; n_dots=0)
    prot[1].name = "Ti"
    prot[1].pdb_element = "Ti"
    @test_throws ArgumentError atomic_sasa(prot)

end
