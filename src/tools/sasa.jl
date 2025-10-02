import Random
import CellListMap
using CellListMap: ParticleSystem, map_pairwise!
using LinearAlgebra: norm
using Statistics: mean

# Container for the custom atom filed that will carry the atom SASA
struct SASA
    sasa::Float32
end

# Structure with the dot cache per atom type
struct DotCache{V}
    x::V
    y::V
    z::V
end

#=
    generate_dots(atomi_radius, probe_radius, n_dots)

Generates points on the surface of a sphere of a given radius (atomic_radius + probe_radius) 
using the double cubic lattice method. `n_dots` controls the maximum number of dots
per sphere.

=#
function generate_dots(atomic_radius, probe_radius::Float32, n_dots::Int; RNG)
    radius = atomic_radius + probe_radius
    if radius <= 0 
        throw(ArgumentError("""\n
            probe_radius too small or incorrectly provided: $probe_radius

        """))
    end
    if n_dots < 1
        throw(ArgumentError("""\n
            n_dots incorrectly set to $n_dots
        """))
    end

    xdot = zeros(Float32, n_dots)
    ydot = zeros(Float32, n_dots)
    zdot = zeros(Float32, n_dots)
    for i in 1:n_dots
        vdot = randn(RNG, SVector{3,Float32})
        vdot = radius * vdot / norm(vdot)
        xdot[i] = vdot[1]
        ydot[i] = vdot[2]
        zdot[i] = vdot[3]
    end

    # Return linear arrays for manual SIMD
    return DotCache(xdot, ydot, zdot)
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

#
# Contribution by Zentrik at:
#
# https://discourse.julialang.org/t/nerd-sniping-can-you-make-this-faster/132793/5?u=lmiq
#
using SIMD: VecRange

function update_dot_exposure!(deltaxy, dot_cache, exposed_i, rj_sq, ::Val{N}) where {N}
    lastN = N * (length(exposed_i) ÷ N)
    lane = VecRange{N}(0)
    @inbounds for i in 1:N:lastN
        if any(exposed_i[lane + i])
            pos_x = dot_cache.x[lane + i] + deltaxy[1]
            pos_y = dot_cache.y[lane + i] + deltaxy[2]
            pos_z = dot_cache.z[lane + i] + deltaxy[3]
            exposed_i[lane + i] &= sum(abs2, (pos_x, pos_y, pos_z)) >= rj_sq
        end
    end
    # Remaining 
    @inbounds for i in lastN+1:length(exposed_i)
        pos_x = dot_cache.x[i] + deltaxy[1]
        pos_y = dot_cache.y[i] + deltaxy[2]
        pos_z = dot_cache.z[i] + deltaxy[3]
        exposed_i[i] &= sum(abs2, (pos_x, pos_y, pos_z)) >= rj_sq
    end
    return exposed_i
end

function update_pair_dot_exposure!(
    x, y, i, j, d2, surface_dots;
    atoms, dot_cache, atom_type::F1, atom_radius_from_type::F2, probe_radius,
    N_SIMD::Val{N}=Val(16),
) where {F1,F2,N}
    type_i = atom_type(atoms[i])
    type_j = atom_type(atoms[j])
    r_i = atom_radius_from_type(type_i) + probe_radius
    r_j = atom_radius_from_type(type_j) + probe_radius
    if d2 <= (r_i + r_j)^2
        R_i_sq = r_i^2
        R_j_sq = r_j^2
        deltaxy = x - y
        update_dot_exposure!(+deltaxy, dot_cache[type_i], surface_dots[i].exposed, R_j_sq, N_SIMD)
        update_dot_exposure!(-deltaxy, dot_cache[type_j], surface_dots[j].exposed, R_i_sq, N_SIMD)
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

- `probe_radius::Real=1.4f0`: The radius of the solvent probe in Angstroms.
- `n_dots::Int=100`: The number of grid points along one axis for dot generation. 
  Higher values lead to more accurate but slower calculations.
- `parallel::Bool=false`: Control if the computation runs in parallel (requires 
  running Julia with multiple threads).

# Example

```julia-repl
julia> using PDBTools

julia> prot = select(read_pdb(PDBTools.TESTPDB), "protein");

julia> at_sasa = atomic_sasa(prot);

julia> sasa(at_sasa) # total sasa of prot
5323.4014f0

julia> sasa(at_sasa, "backbone") # backbone sasa in prot
983.07745f0

julia> sasa(at_sasa, "not backbone") # other atoms
4442.001f0

julia> sasa(at_sasa, "resname ARG GLU") # some residue types
534.0886f0
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
    probe_radius::Real=1.4f0,
    n_dots::Int=100,
    atom_type::Function=element,
    atom_radius_from_type::Function=type -> getproperty(elements[type], :vdw_radius),
    parallel=false,
    N_SIMD::Val{N}=Val(16), # Size of SIMD blocks. Can be tunned for maximum performance.
    RNG=Random.default_rng(), # used for testing with StableRNGs
) where {N}
    probe_radius = Float32(probe_radius)

    # Unique list of atom types
    atom_types = atom_type.(unique(atom_type, atoms))

    # Memoization for dot generation to avoid recomputing for same radii
    dot_cache = Dict{eltype(atom_types),DotCache{Vector{Float32}}}()
    for type in atom_types
        atom_radius = atom_radius_from_type(type)
        if isnan(atom_radius)
            throw(ArgumentError("""\n
                Atom of type $type does not have a vdW radius defined.
                Use custom `atom_type` and `atom_radius_from_type` input parameters if needed. 
            """))
        end
        dot_cache[type] = generate_dots(atom_radius_from_type(type), probe_radius, n_dots; RNG)
    end

    system = ParticleSystem(
        xpositions=coor.(atoms),
        unitcell=nothing,
        cutoff=2 * (maximum(atom_radius_from_type(type) for type in atom_types) + probe_radius),
        output=[AtomDots(trues(length(dot_cache[atom_type(at)].x))) for at in atoms],
        output_name=:surface_dots,
        parallel=parallel,
    )

    map_pairwise!(
        (x, y, i, j, d2, surface_dots) ->
            update_pair_dot_exposure!(
                x, y, i, j, d2, surface_dots;
                atoms, dot_cache, atom_type, atom_radius_from_type, probe_radius, N_SIMD,
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

```julia-repl
julia> using PDBTools

julia> prot = select(read_pdb(PDBTools.TESTPDB), "protein");

julia> at_sasa = atomic_sasa(prot);

julia> sasa(at_sasa) # total sasa of prot
5325.885f0

julia> sasa(at_sasa, "backbone") # selection string
971.95325f0

julia> sasa(at_sasa, at -> name(at) == "CA") # selection function
45.88987f0

julia> sasa(at_sasa[1]) # single atom SASA
5.467941f0
```

"""
function sasa end

sasa(atom::Atom{SASA}) = atom.custom.sasa
sasa(atoms::AbstractVector{<:Atom{SASA}}) = sum(sasa, atoms)
sasa(atoms::AbstractVector{<:Atom{SASA}}, sel::String) = sasa(atoms, Select(sel))
sasa(atoms::AbstractVector{<:Atom{SASA}}, sel::Function) = sum(sasa(at) for at in atoms if sel(at))

@testitem "sasa" begin
    using StableRNGs
    using PDBTools
    prot = read_pdb(PDBTools.TESTPDB, "protein")
    N = 10_000 # number of dot samples
    RNG=StableRNG(123)

    # Compare with the output of VMD: the difference is that VMD uses a vdW radius for H 
    # of 1.00, and we use 1.10.
    # vmd version 2.0.0a5 was run with e. g.
    # [ measure sasa 1.4 $protein -samples 100000 -restrict [ atomselect top "backbone" ] ]
    const vmd_radii = Dict(
        "N" => 1.55,
        "C" => 1.70,
        "H" => 1.00,
        "O" => 1.52,
        "S" => 1.80,
    )
    at_sasa = atomic_sasa(prot; n_dots=N, atom_radius_from_type=type -> vmd_radii[type], RNG)
    @test sasa(at_sasa) ≈ 5365.55029296875 rtol=0.01
    # Accessiblity of groups within the structure
    @test sasa(at_sasa, "backbone") ≈ 1130.37646484375 rtol=0.05
    @test sasa(at_sasa, "resname GLU LYS") ≈ 797.8261108398438 rtol=0.05
    @test sasa(at_sasa, "residue 1") ≈ 124.57905578613281 rtol=0.05
    @test sasa(at_sasa, "residue 104") ≈ 122.50507354736328 rtol=0.05

    # Compare with Gromacs - 2023.3 output
    # gmx sasa -s prot.pdb -o sasa_output.xvg -ndots 100000
    @test sasa(atomic_sasa(prot; n_dots=N, RNG)) ≈ 100 * 53.754 rtol=0.01
    # Isolated groups
    @test sasa(atomic_sasa(select(prot, "backbone and not name O"); n_dots=N, RNG)) ≈ 100 * 55.229 rtol = 0.05
    @test sasa(atomic_sasa(select(prot, "name CA"); n_dots=N, RNG)) ≈ 100 * 58.630 rtol = 0.05
    @test sasa(atomic_sasa(select(prot, "sidechain and not element H"); n_dots=N, RNG)) ≈ 100 * 69.024 rtol = 0.05

    # Test non-contiguous indexing with general selections
    at_sasa = atomic_sasa(select(prot, "name CA"); n_dots=N, RNG)
    @test sasa(at_sasa) ≈ 5860.0 rtol=0.05
    @test sasa(at_sasa, "resname THR") ≈ 322.2 rtol=0.05
    @test sasa(at_sasa, at -> resname(at) == "THR") ≈ 322.2 rtol=0.05

    # Test parallelization
    @test sasa(atomic_sasa(prot; parallel=false, n_dots=N)) ≈ sasa(atomic_sasa(prot; parallel=true, n_dots=N)) rtol=0.05

    # Test errors
    @test_throws ArgumentError atomic_sasa(prot; probe_radius=-2.0)
    @test_throws ArgumentError atomic_sasa(prot; n_dots=0)
    prot[1].name = "Ti"
    prot[1].pdb_element = "Ti"
    @test_throws ArgumentError atomic_sasa(prot)

end
