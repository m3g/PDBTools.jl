using LinearAlgebra

"""
Simple structure to hold atomic coordinates and radius.
"""
struct Atom
    x::Float64
    y::Float64
    z::Float64
    radius::Float64
end

"""
    generate_dots(radius, n_grid_points)

Generates points on the surface of a sphere of a given `radius` using the 
double cubic lattice method. `n_grid_points` controls the density of the points.
"""
function generate_dots(radius::Float64, n_grid_points::Int)
    dots = Vector{Vector{Float64}}()
    if radius <= 0.0 || n_grid_points <= 0
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
            push!(dots, [x, y, z] .* norm_factor)
        end
    end

    # --- Second Lattice (shifted by half a step) ---
    offset = step / 2.0
    for i in 0:n_grid_points, j in 0:n_grid_points, k in 0:n_grid_points
        x = -radius + offset + i * step
        y = -radius + offset + j * step
        z = -radius + offset + k * step
        dist_sq = x^2 + y^2 + z^2
        
        # Project internal points onto the surface
        if dist_sq < radius_sq
            norm_factor = radius / sqrt(dist_sq)
            push!(dots, [x, y, z] .* norm_factor)
        end
    end
    
    # Return unique points to avoid redundancy
    unique_dots_set = Set(Tuple(round.(d, digits=8)) for d in dots)
    return [collect(t) for t in unique_dots_set]
end


"""
    calculate_sasa(atoms; probe_radius, n_grid_points)

Calculates the Solvent Accessible Surface Area (SASA) for a vector of `Atom`s.

# Arguments
- `atoms::Vector{Atom}`: A vector of atoms in the molecule.
- `probe_radius::Float64=1.4`: The radius of the solvent probe in Angstroms.
- `n_grid_points::Int=20`: The number of grid points along one axis for dot generation. 
  Higher values lead to more accurate but slower calculations.
"""
function calculate_sasa(atoms::Vector{Atom}; probe_radius::Float64 = 1.4, n_grid_points::Int = 20)
    total_sasa = 0.0
    
    # Create atoms with their radii extended by the probe radius
    extended_atoms = [Atom(a.x, a.y, a.z, a.radius + probe_radius) for a in atoms]
    
    # Memoization for dot generation to avoid recomputing for same radii
    dot_cache = Dict{Float64, Vector{Vector{Float64}}}()

    for (i, atom_i) in enumerate(extended_atoms)
        R_i = atom_i.radius
        center_i = [atom_i.x, atom_i.y, atom_i.z]
        
        # Generate or retrieve the dots for this atom's extended radius
        if !haskey(dot_cache, R_i)
            dot_cache[R_i] = generate_dots(R_i, n_grid_points)
        end
        dots_template = dot_cache[R_i]
        
        n_total_dots = length(dots_template)
        if n_total_dots == 0
            continue
        end
        
        n_accessible_dots = 0
        for dot in dots_template
            # Position the dot on the atom's surface in the molecule's coordinate system
            dot_on_surface = center_i + dot
            
            is_occluded = false
            for (j, atom_j) in enumerate(extended_atoms)
                if i == j
                    continue
                end
                
                center_j = [atom_j.x, atom_j.y, atom_j.z]
                R_j = atom_j.radius
                
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
        atom_sasa = (n_accessible_dots / n_total_dots) * 4.0 * π * R_i^2
        total_sasa += atom_sasa
    end
    
    return total_sasa
end

# --- Example Usage ---
# Define two carbon atoms (vdW radius ≈ 1.7 Å) separated by 3.0 Å
atom1 = Atom(0.0, 0.0, 0.0, 1.7)
atom2 = Atom(3.0, 0.0, 0.0, 1.7)
molecule = [atom1, atom2]

# Calculate SASA
sasa_value = calculate_sasa(molecule, probe_radius=1.4, n_grid_points=25)

# For comparison, the SASA of a single, isolated atom
isolated_atom_sasa = calculate_sasa([atom1], probe_radius=1.4, n_grid_points=25)


println("SASA of the two-atom molecule: ", round(sasa_value, digits=2), " Å²")
println("SASA of a single isolated atom: ", round(isolated_atom_sasa, digits=2), " Å²")
println("Expected sum for two isolated atoms: ", round(2 * isolated_atom_sasa, digits=2), " Å²")