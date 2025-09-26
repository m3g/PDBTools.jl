using LinearAlgebra

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


"""
    sasa(atoms; probe_radius, n_grid_points)

Calculates the Solvent Accessible Surface Area (SASA) for a vector of `Atom`s.

# Arguments
- `atoms::Vector{Atom}`: A vector of atoms in the molecule.
- `probe_radius::Float64=1.4`: The radius of the solvent probe in Angstroms.
- `n_grid_points::Int=20`: The number of grid points along one axis for dot generation. 
  Higher values lead to more accurate but slower calculations.
"""
function sasa(atoms::AbstractVector{<:Atom}; probe_radius::Real = 1.4, n_grid_points::Int = 20)
    total_sasa = 0.0f0
    
    # Memoization for dot generation to avoid recomputing for same radii
    dot_cache = Dict{String, Vector{Vector{Float32}}}()

    for (i, atom_i) in enumerate(atoms)
        el_i = element(atom_i)
        R_i = element_vdw_radius(atom_i) + probe_radius
        center_i = SVector{3,Float32}(atom_i.x, atom_i.y, atom_i.z)
        
        # Generate or retrieve the dots for this atom's extended radius
        if !haskey(dot_cache, el_i)
            dot_cache[el_i] = generate_dots(R_i, n_grid_points)
        end
        dots_template = dot_cache[el_i]
        
        n_total_dots = length(dots_template)
        if n_total_dots == 0
            continue
        end
        
        n_accessible_dots = 0
        for dot in dots_template
            # Position the dot on the atom's surface in the molecule's coordinate system
            dot_on_surface = center_i + dot
            
            is_occluded = false
            for (j, atom_j) in enumerate(atoms)
                if i == j
                    continue
                end
                
                center_j = SVector{3,Float32}(atom_j.x, atom_j.y, atom_j.z)
                R_j = element_vdw_radius(atom_j)
                
                # Check if the dot is inside the neighboring atom j
                if sum((dot_on_surface - center_j).^2) < R_j^2
                    is_occluded = true
                    break # Dot is buried, move to the next dot
                end
            end
            
            if !is_occluded
                n_accessible_dots += 1
            end
        end
        
        # Calculate and add this atom's contribution to the total SASA
        atom_sasa = 4 * π * R_i^2 * (n_accessible_dots / n_total_dots) 
        total_sasa += atom_sasa
    end
    
    return total_sasa
end

@testitem "sasa" begin

    sasa_value = calculate_sasa(molecule, probe_radius=1.4, n_grid_points=25)

    # For comparison, the SASA of a single, isolated atom
    isolated_atom_sasa = calculate_sasa([atom1], probe_radius=1.4, n_grid_points=25)
    
    
    println("SASA of the two-atom molecule: ", round(sasa_value, digits=2), " Å²")
    println("SASA of a single isolated atom: ", round(isolated_atom_sasa, digits=2), " Å²")
    println("Expected sum for two isolated atoms: ", round(2 * isolated_atom_sasa, digits=2), " Å²")

end