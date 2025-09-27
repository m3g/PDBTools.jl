using LinearAlgebra

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

function update_dot_oclusion!(i, j, x, y, atoms, dot_cache, system)
    el_i = element(atoms[i])
    R_j_sq = element_vdw_radius(atoms[j])^2
    for idot in enumerate(system.ocluded_atom_dots[i])
        if !ocluded_atom_dots[i]
            # Position the dot on the atom's surface in the molecule's coordinate system
            dot_on_surface = x + dot_cache[el_i][idot]
            # Check if the dot is inside the neighboring atom j
            if sum((dot_on_surface - y).^2) < R_j_sq
                system.ocluded_atom_dots[i] = true
            end
        end
    end
    return nothing
end

function pair_dot_oclusion!(x, y, i, j, ocluded_atom_dots; atoms, dot_cache, system)
    update_dot_oclusion!(i, j, x, y, atoms, dot_cache, system)
    update_dot_oculsion!(j, i, y, x, atoms, dot_cache, system)
    return ocluded_atom_dots
end

"""
    sasa(atoms; probe_radius, n_grid_points)

Calculates the Solvent Accessible Surface Area (SASA) for a vector of `Atom`s.

# Arguments
- `atoms::Vector{Atom}`: A vector of atoms in the molecule.
- `probe_radius::Float64=1.4`: The radius of the solvent probe in Angstroms.
- `n_grid_points::Int=20`: The number of grid points along one axis for dot generation. 
  Higher values lead to more accurate but slower calculations.
"""
function compute_sasa(atoms::AbstractVector{<:Atom}; probe_radius::Real = 1.4, n_grid_points::Int = 20)

    elms = unique(element, atoms)
    radius = [ elements[el].vdw_radius for el in elms ]
    
    # Memoization for dot generation to avoid recomputing for same radii
    dot_cache = Dict{String, Vector{Vector{Float32}}}()
    for (iel, el) in enumerate(elms)
        dot_cache[el] = generate_dots(radius[iel], n_grid_points)
    end

    sys = ParticleSystem(
        xpositions = coor.(atoms),
        unitcell = nothing,
        cutoff = maximum(radius) + probe_radius,
        output = [ zeros(Bool, length(dot_cache[element(at)])) for at in atoms ],
        output_name = :ocluded_atom_dots
        parallel=false,
    )

    map_pairwise(
        (x, y, i, j, d2, ocluded_atom_dots) -> 
            pair_dot_oclusion!(x, y, i, j, ocluded_atom_dots; atoms, dot_cache, system),
        sys
    )

    ats_with_sasa = [
        add_custom_field(
            atoms[i], 
            SASA(4*π*element_vdw_radius(atoms[i])^2*(count(==(false), system.ocluded_atom_dots)/n_grid_points))
        )
    ]
    
    return ats_with_sasa
end

function sasa(atoms::AbstractVector{<:Atom{SASA}})
    return sum(at.custom.sasa for at in atoms)
end

function sasa(atoms::AbstractVector{<:Atom{SASA}}, selection::Union{Function,String})
    isel = index.(select(atoms, selection))
    return sum(atoms[i].custom.sasa for i in isel)
end

@testitem "sasa" begin

    sasa_value = calculate_sasa(molecule, probe_radius=1.4, n_grid_points=25)

    # For comparison, the SASA of a single, isolated atom
    isolated_atom_sasa = calculate_sasa([atom1], probe_radius=1.4, n_grid_points=25)
    
    
    println("SASA of the two-atom molecule: ", round(sasa_value, digits=2), " Å²")
    println("SASA of a single isolated atom: ", round(isolated_atom_sasa, digits=2), " Å²")
    println("Expected sum for two isolated atoms: ", round(2 * isolated_atom_sasa, digits=2), " Å²")

end
