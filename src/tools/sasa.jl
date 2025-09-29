import CellListMap
using CellListMap: ParticleSystem, map_pairwise!
using LinearAlgebra
using Statistics: mean

struct SASA sasa::Float32 end

"""
    generate_dots(radius, n_grid_points)

Generates points on the surface of a sphere of a given `radius` using the 
double cubic lattice method. `n_grid_points` controls the density of the points.
"""
function generate_dots(radius::Real, n_grid_points::Int)
    dots = SVector{3,Float32}[]
    if radius <= 0 || n_grid_points <= 0
        return dots
    end
    
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

struct AtomDots0
    exposed::Vector{Bool}
end
function CellListMap.reset_output!(x::AtomDots0) 
    x.exposed .= true
    return x
end
CellListMap.copy_output(x::AtomDots0) = AtomDots0(copy(x.exposed))
function CellListMap.reducer(x::AtomDots0, y::AtomDots0) 
    x.exposed .= x.exposed .& y.exposed
    return x
end

function update_dot_oclusion!(
    i, j, x, y, atoms, dot_cache, surface_dots; 
    atom_type::Function,
    atom_radius_from_type::Function,
    probe_radius,
)
    type_i = atom_type(atoms[i])
    R_j_sq = (atom_radius_from_type(atom_type(atoms[j])) + probe_radius)^2
    for (idot, dot_exposed) in enumerate(surface_dots[i].exposed)
        if dot_exposed
            # Position the dot on the atom's surface in the molecule's coordinate system
            dot_on_surface = x + dot_cache[type_i][idot]
            # Check if the dot is inside the neighboring atom j
            if sum(abs2, dot_on_surface - y) < R_j_sq
                surface_dots[i].exposed[idot] = false
            end
        end
    end
    return nothing
end

function pair_dot_oclusion!(
    x, y, i, j, surface_dots; 
    atoms, dot_cache, atom_type::F1, atom_radius_from_type::F2, probe_radius,
) where {F1,F2}
    update_dot_oclusion!(i, j, x, y, atoms, dot_cache, surface_dots; atom_type, atom_radius_from_type, probe_radius)
    update_dot_oclusion!(j, i, y, x, atoms, dot_cache, surface_dots; atom_type, atom_radius_from_type, probe_radius)
    return surface_dots
end

"""
    sasa(atoms; probe_radius, n_grid_points)

Calculates the Solvent Accessible Surface Area (SASA) for a vector of `Atom`s.

# Arguments

- `atoms::Vector{PDBTools.Atom}`: A vector of atoms in the molecule.
- `probe_radius::Real=1.4`: The radius of the solvent probe in Angstroms.
- `n_grid_points::Int=20`: The number of grid points along one axis for dot generation. 
  Higher values lead to more accurate but slower calculations.

"""
function compute_sasa(
    atoms::AbstractVector{<:Atom}; 
    probe_radius::Real = 1.4, 
    n_grid_points::Int = 20,
    atom_type::Function = element,
    atom_radius_from_type::Function = type -> getproperty(elements[type], :vdw_radius),
    parallel=true,
)

    # Unique list of atom types
    atom_types = atom_type.(unique(atom_type, atoms))
    
    # Memoization for dot generation to avoid recomputing for same radii
    dot_cache = Dict{eltype(atom_types), Vector{SVector{3,Float32}}}()
    for type in atom_types
        dot_cache[type] = generate_dots(atom_radius_from_type(type) + probe_radius, n_grid_points)
    end
        
    system = ParticleSystem(
        xpositions = coor.(atoms),
        unitcell = nothing,
        cutoff = 2*maximum(atom_radius_from_type(type) for type in atom_types) + probe_radius,
        output = [ AtomDots0(ones(Bool, length(dot_cache[atom_type(at)]))) for at in atoms ],
        output_name = :surface_dots,
        parallel=parallel,
    )

    map_pairwise!(
        (x, y, i, j, d2, surface_dots) -> 
            pair_dot_oclusion!(
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

export compute_sasa, sasa

sasa(atom::Atom{SASA}) = atom.custom.sasa
sasa(atoms::AbstractVector{<:Atom{SASA}}) = sum(sasa, atoms)
function sasa(atoms::AbstractVector{<:Atom{SASA}}, selection::Union{Function,String})
    isel = index.(select(atoms, selection))
    return sum(sasa(atoms[i]) for i in isel)
end

@testitem "sasa" begin

    # gmx sasa -s prot.pdb -probe 1.4 -ndots 20 -n index.ndx -surface CYSTHR -output CYSTHR
    # Using the apparent vdW radii that Gromacs uses
    gmx_radii = Dict(
        "N" => 1.555,
        "C" => 1.570,
        "H" => 1.520,
        "O" => 1.552,
        "S" => 1.575,
    )
    atom_radius_from_type(type) = gmx_radii[type]

    # gmx sasa -s prot.pdb -probe 1.4 -ndots 20 -n index.ndx -surface CYSTHR -output CYSTHR
    # Using the apparent vdW radii that Gromacs uses
    vmd_radii = Dict(
        "N" => 1.55,
        "C" => 1.70,
        "H" => 1.00,
        "O" => 1.52,
        "S" => 1.80,
    )
    atom_radius_from_type(type) = vmd_radii[type]

    sasa_gmx = Dict(
        "total" => 115.384,
        "first_atom" => 0.0,
        "1st_and_last (isolated)" => 43.657,
        "first_residue (isolated)" => 34.182,
        "CYS and THR" => 7.734,
        "CYS and THR (isolated)" => 80.959,
        "N atom (isolated)" => 30.386,
        "S atom (isolated)" => 31.171,
        "C atom (isolated)" => 30.975,
        "O atom (isolated)" => 30.269,
        "H atom (isolated)" => 29.033,
    )

    prot = read_pdb(PDBTools.TESTPDB, "protein")
    @test sasa(compute_sasa([prot[1]])) ≈ sasa_gmx["first atom (isolated)"] rtol = 1e-3
    @test sasa(compute_sasa([prot[20]])) ≈ sasa_gmx["S atom (isolated)"] rtol = 1e-3
    @test sasa(compute_sasa([prot[5]])) ≈ sasa_gmx["C atom (isolated)"] rtol = 1e-3
    @test sasa(compute_sasa([prot[12]])) ≈ sasa_gmx["C atom (isolated)"] rtol = 1e-3
    @test sasa(compute_sasa([prot[2]])) ≈ sasa_gmx["H atom (isolated)"] rtol = 1e-3


    # For comparison, the SASA of a single, isolated atom
    isolated_atom_sasa = calculate_sasa([atom1], probe_radius=1.4, n_grid_points=25)
    
    
    println("SASA of the two-atom molecule: ", round(sasa_value, digits=2), " Å²")
    println("SASA of a single isolated atom: ", round(isolated_atom_sasa, digits=2), " Å²")
    println("Expected sum for two isolated atoms: ", round(2 * isolated_atom_sasa, digits=2), " Å²")

end
