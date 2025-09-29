import CellListMap
using CellListMap: ParticleSystem, map_pairwise!
using LinearAlgebra
using Statistics: mean

# Container for the custom atom filed that will carry the atom SASA
struct SASA sasa::Float32 end

#=
    generate_dots(atomi_radius, probe_radius, n_dots)

Generates points on the surface of a sphere of a given radius (atomic_radius + probe_radius) 
using the double cubic lattice method. `n_dots` controls the maximum number of dots
per sphere.

=#
function generate_dots(atomic_radius, probe_radius::Real, n_dots::Int)
    radius = atomic_radius + probe_radius
    n_grid_points = round(Int, (3/(4 * π * ((1/2)^3)) * n_dots)^(1/3))
    if radius <= 0 || n_grid_points <= 0
        throw(ArgumentError("""\n
            The probe_radius ($probe_radius) or number of dots ($n_dots) are too small, or incorrectly provided.

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
    exposed::Vector{Bool}
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

function update_dot_exposure!(
    i, j, x, y, atoms, dot_cache, surface_dots; 
    atom_type::Function,
    atom_radius_from_type::Function,
    probe_radius,
)
    type_i = atom_type(atoms[i])
    type_j = atom_type(atoms[j])
    R_j_sq = (atom_radius_from_type(type_j) + probe_radius)^2
    for (idot, dot_exposed) in enumerate(surface_dots[i].exposed)
        if dot_exposed
            dot_on_surface = x + dot_cache[type_i][idot]
            # Position the dot on the atom's surface in the molecule's coordinate system
            # Check if the dot is inside the neighboring atom j
            if sum(abs2, dot_on_surface - y) < R_j_sq
                surface_dots[i].exposed[idot] = false
            end
        end
    end
    return nothing
end

function update_pair_dot_exposure!(
    x, y, i, j, surface_dots; 
    atoms, dot_cache, atom_type::F1, atom_radius_from_type::F2, probe_radius,
) where {F1,F2}
    update_dot_exposure!(i, j, x, y, atoms, dot_cache, surface_dots; atom_type, atom_radius_from_type, probe_radius)
    update_dot_exposure!(j, i, y, x, atoms, dot_cache, surface_dots; atom_type, atom_radius_from_type, probe_radius)
    return surface_dots
end

"""
    atomic_sasa(atoms; probe_radius, n_dots)::Vector{PDBTools.Atom{SASA}}

Calculates the Solvent Accessible Surface Area (SASA) for a vector of `Atom`s. 
Returns a vector of `PDBToools.Atom`s, with a custom fiels with `PDBTools.SASA` 
structure, which contains the solvent accessible surface area for that atom.

The `sasa` function computes the total SASA or the SASA of a subset of the atoms
in the structure. The output value is provided in Å^2. 

# Main argument

- `atoms::Vector{PDBTools.Atom}`: A vector of atoms in the molecule.

# Optional arguments 

- `probe_radius::Real=1.4`: The radius of the solvent probe in Angstroms.
- `n_dots::Int=500`: The number of grid points along one axis for dot generation. 
  Higher values lead to more accurate but slower calculations.
- `parallel::Bool=true`: Control if the computation runs in parallel (requires 
  running Julia with multiple threads).

# Example

```jldoctest
julia> using PDBTools

julia> prot = select(read_pdb(PDBTools.TESTPDB), "protein");

julia> at_sasa = atomic_sasa(prot);

julia> sasa(at_sasa) # total sasa of prot
5374.4414f0

julia> sasa(at_sasa, "backbone") # backbone sasa in prot
873.72626f0

julia> sasa(at_sasa, "not backbone") # other atoms
4500.7124f0

julia> sasa(at_sasa, "resname ARG GLU") # some residue types
531.51776f0
```

# Additional control:

Two arguments can control the atom radii used for computing the SASA. These arguments
are functions:

- `atom_type`: Function that given each atom of the array of atoms, returns the atom "type".
- `atom_radius_from_type`: Given the atom "type", returns the vdW radius of the atom. 

By default, `atom_type = PDBTools.element`, a function that just returns the element of the atom,
and `atom_radius_from_type` obtains the vdW radius from the `PDBTools.elements` list, in which
the atomc radii of https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page) are See the documentation for further information.
used. Atoms with missing radius have a `NaN` value, and the computation will not return meaningful
values. See the documentation for further information.

"""
function atomic_sasa(
    atoms::AbstractVector{<:Atom}; 
    probe_radius::Real=1.4, 
    n_dots::Int=500,
    atom_type::Function = element,
    atom_radius_from_type::Function = type -> getproperty(elements[type], :vdw_radius),
    parallel=true,
)

    # Unique list of atom types
    atom_types = atom_type.(unique(atom_type, atoms))
    
    # Memoization for dot generation to avoid recomputing for same radii
    dot_cache = Dict{eltype(atom_types), Vector{SVector{3,Float32}}}()
    for type in atom_types
        dot_cache[type] = generate_dots(atom_radius_from_type(type), probe_radius, n_dots)
    end
        
    system = ParticleSystem(
        xpositions = coor.(atoms),
        unitcell = nothing,
        cutoff = 2*(maximum(atom_radius_from_type(type) for type in atom_types) + probe_radius),
        output = [ AtomDots(ones(Bool, length(dot_cache[atom_type(at)]))) for at in atoms ],
        output_name = :surface_dots,
        parallel=parallel,
    )

    map_pairwise!(
        (x, y, i, j, d2, surface_dots) -> 
            update_pair_dot_exposure!(
                x, y, i, j, surface_dots; 
                atoms, dot_cache, atom_type, atom_radius_from_type, probe_radius,
            ),
        system,
    )

    ats_with_sasa = [
        add_custom_field(
            atoms[i], 
            SASA( 
                4*π*(atom_radius_from_type(atom_type(atoms[i])) + probe_radius)^2 *
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

# Example

```jldoctest
julia> using PDBTools

julia> prot = select(read_pdb(PDBTools.TESTPDB), "protein");

julia> at_sasa = atomic_sasa(prot);

julia> sasa(at_sasa) # total sasa of prot
5374.4414f0

julia> sasa(at_sasa, "backbone") # selection string
873.72626f0

julia> sasa(at_sasa, at -> name(at) == "CA") # selection function
20.5797f0

julia> sasa(at_sasa[1]) # single atom SASA
1.8389939f0
```

"""
function sasa end

sasa(atom::Atom{SASA}) = atom.custom.sasa
sasa(atoms::AbstractVector{<:Atom{SASA}}) = sum(sasa, atoms)
function sasa(atoms::AbstractVector{<:Atom{SASA}}, selection::Union{Function,String})
    isel = index.(select(atoms, selection))
    return sum(sasa(atoms[i]) for i in isel)
end

@testitem "sasa" begin

    # Using VMD radii
    vmd_radii = Dict(
        "N" => 1.55,
        "C" => 1.70,
        "H" => 1.00,
        "O" => 1.52,
        "S" => 1.80,
    )

    prot = read_pdb(PDBTools.TESTPDB, "protein")

    @test sasa(atomic_sasa([prot[1]])) ≈ sasa_gmx["first atom (isolated)"] rtol = 1e-3
    @test sasa(atomic_sasa([prot[20]])) ≈ sasa_gmx["S atom (isolated)"] rtol = 1e-3
    @test sasa(atomic_sasa([prot[5]])) ≈ sasa_gmx["C atom (isolated)"] rtol = 1e-3
    @test sasa(atomic_sasa([prot[12]])) ≈ sasa_gmx["C atom (isolated)"] rtol = 1e-3
    @test sasa(atomic_sasa([prot[2]])) ≈ sasa_gmx["H atom (isolated)"] rtol = 1e-3


    # For comparison, the SASA of a single, isolated atom
    isolated_atom_sasa = atomic_sasa([atom1], probe_radius=1.4)
    
    
    println("SASA of the two-atom molecule: ", round(sasa_value, digits=2), " Å²")
    println("SASA of a single isolated atom: ", round(isolated_atom_sasa, digits=2), " Å²")
    println("Expected sum for two isolated atoms: ", round(2 * isolated_atom_sasa, digits=2), " Å²")

end
