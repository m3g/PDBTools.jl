import CellListMap
using CellListMap: ParticleSystem, pairwise!
using LinearAlgebra: norm
using SIMD: VecRange

abstract type AtomicRadiiType end
struct CustomAtomicRadii <: AtomicRadiiType end
struct StandardAtomicRadii <: AtomicRadiiType end

#
# Container for the custom atom filed that will carry the resulting SASA information
#
struct SASA{RadiiType<:AtomicRadiiType,N,V}
    particles::V
    sasa::Vector{Float32}
    dots::Vector{Vector{SVector{N,Float32}}}
end
Base.getindex(x::SASA, i::Integer) = x.sasa[i]
Base.eachindex(x::SASA) = eachindex(x.sasa)
function Base.show(io::IO, ::MIME"text/plain", s::SASA)
    print(io, chomp("""
    $(typeof(s))
        Number of particles: $(length(s.particles))
        Total SASA: $(sasa(s))
        Output of dots: $(isempty(s.dots) ? false : true) 
    """))
end

# Structure with the dot cache per atom type
struct DotCache{V}
    x::V
    y::V
    z::V
end

#=
    generate_dots(atomic_radius::T, probe_radius::T, n_dots::Int) where {T<:Real}

Generates points on the surface of a sphere of a given radius (atomic_radius + probe_radius) 
using a Fibonacci lattice. The type of `probe_radius` defines the output type of the 
coordinates.

=#
function generate_dots(atomic_radius::Real, probe_radius::T, n_dots::Int) where {T<:Real}
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
    # Return linear arrays for manual SIMD
    xdot = zeros(T, n_dots)
    ydot = zeros(T, n_dots)
    zdot = zeros(T, n_dots)
    phi = (1 + sqrt(5)) / 2
    for i in 1:n_dots
        theta = acos(1 - 2 * i / n_dots)
        phi_angle = 2 * pi * i / phi
        xdot[i] = radius * cos(phi_angle) * sin(theta)
        ydot[i] = radius * sin(phi_angle) * sin(theta)
        zdot[i] = radius * cos(theta)
    end
    return DotCache(xdot, ydot, zdot)
end

#
# Structure that carries the information about dots of each atom, if 
# they are exposed or found to be occluded by other atoms.
#
struct AtomDots{T<:AbstractVector{Bool}}
    exposed::T
end
# Make a contiguous chunk of memory for all exposed information
struct AtomDotMatrix
    all_dots::Matrix{Bool}
end
Base.getindex(m::AtomDotMatrix, i) = AtomDots(@view(m.all_dots[:, i]))
CellListMap.copy_output(m::AtomDotMatrix) = AtomDotMatrix(copy(m.all_dots))
function CellListMap.reset_output!(m::AtomDotMatrix)
    m.all_dots .= true
    return m
end
function CellListMap.reducer(x::AtomDotMatrix, y::AtomDotMatrix)
    x.all_dots .= x.all_dots .& y.all_dots
    return x
end

#
# Contribution by Zentrik at:
#
# https://discourse.julialang.org/t/nerd-sniping-can-you-make-this-faster/132793/5?u=lmiq
#
function update_dot_exposure!(deltaxy, dot_cache, exposed_i, rj_sq, ::Val{N}) where {N}
    lastN = N * (length(exposed_i) ÷ N)
    lane = VecRange{N}(0)
    @inbounds for i in 1:N:lastN
        if any(exposed_i[lane+i])
            pos_x = dot_cache.x[lane+i] + deltaxy[1]
            pos_y = dot_cache.y[lane+i] + deltaxy[2]
            pos_z = dot_cache.z[lane+i] + deltaxy[3]
            exposed_i[lane+i] &= sum(abs2, (pos_x, pos_y, pos_z)) >= rj_sq
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
    pair, surface_dots;
    atoms, dot_cache, atom_type::F1, atom_radius_from_type::F2, probe_radius,
    N_SIMD::Val{N}=Val(16),
) where {F1,F2,N}
    (; x, y, i, j, d2) = pair
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
    sasa_particles(atoms; probe_radius, n_dots)

Calculates the Solvent Accessible Surface Area (SASA) for a vector of `Atom`s. 

# Main argument

- `atoms::Vector{PDBTools.Atom}`: A vector of atoms in the molecule.

# Returns

`PDBTools.SASA` structure, containing the vector of atoms, the SASA of each atom (in Å²) and,
optionally, the solvent accessible dots that define the surface. 

The `sasa` function computes the total SASA or the SASA of a subset of the atoms
in the structure.

# Optional arguments 

- `probe_radius::Real=1.4f0`: The radius of the solvent probe in Angstroms.
- `n_dots::Int=512`: The number of grid points along one axis for dot generation. 
  Higher values lead to more accurate but slower calculations.
- `unitcell=nothing`: if periodic boundary conditions are used, provide a 3x3 matrix with
  the unitcell, or alternatively a vector of length 3 with the sides, for orthorhombic cells.
- `parallel::Bool=true`: Control if the computation runs in parallel (requires 
  running Julia with multiple threads).

# Example

```julia-repl
julia> using PDBTools

julia> prot = select(read_pdb(PDBTools.TESTPDB), "protein");

julia> at_sasa = sasa_particles(prot);

julia> sasa(at_sasa) # total sasa of prot
5389.0146f0

julia> sasa(at_sasa, "backbone") # backbone sasa in prot
988.7648f0

julia> sasa(at_sasa, "not backbone") # other atoms
4400.246f0

julia> sasa(at_sasa, "resname ARG GLU") # some residue types
543.29846f0
```

# Additional control:

Two arguments can control the atom radii used for computing the SASA. These arguments
are functions:

- `atom_type`: Function that given each atom of the array of atoms, returns the atom "type".
- `atom_radius_from_type`: Given the atom "type", returns the vdW radius of the atom. 
- `output_dots::Bool=false`: If true, the resulting `SASA` structure will contain the solvent accessible
    dots per particle in the `dots` field.

By default, `atom_type = PDBTools.element`, a function that just returns the element symbol of the atom,
and `atom_radius_from_type` obtains the vdW radius from the `PDBTools.elements` list given the element symbol.
Here, the atomic radii of https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page) are used. 
Atoms with missing radius have a `NaN` value, and the computation will not return meaningful
values. 

"""
function sasa_particles(
    atoms::AbstractVector{<:Atom};
    probe_radius::Real=1.4f0,
    n_dots::Int=512,
    atom_type::Union{Nothing,Function}=nothing,
    atom_radius_from_type::Union{Nothing,Function}=nothing,
    output_dots::Bool=false,
    unitcell::Union{AbstractVector,AbstractMatrix,Nothing}=nothing,
    parallel=true,
    N_SIMD::Val{N}=Val(16), # Size of SIMD blocks. Can be tuned for maximum performance.
) where {N}
    # Without defining atom type functions, return default StandardAtomicRadii calculation
    if isnothing(atom_type) & isnothing(atom_radius_from_type)
        return sasa_particles(StandardAtomicRadii, atoms;
            probe_radius, n_dots, output_dots, unitcell, parallel, N_SIMD
        )
    else
        # Return a CustomAtomicRadii structure
        atom_type = isnothing(atom_type) ? element : atom_type
        atom_radius_from_type = isnothing(atom_radius_from_type) ? 
            type -> getproperty(elements[type], :vdw_radius) : atom_radius_from_type
        return _sasa_particles(CustomAtomicRadii, atoms;
            atom_type=atom_type,
            atom_radius_from_type=atom_radius_from_type,
            probe_radius, n_dots, output_dots, unitcell, parallel, N_SIMD
        )
    end
end 

#
# Internal main sasa calculation interface 
#
function _sasa_particles(
    RadiiType::Type{<:AtomicRadiiType},
    atoms::AbstractVector{<:Atom};
    probe_radius::Real=1.4f0,
    n_dots::Int=512,
    atom_type::Function=element,
    atom_radius_from_type::Function=type -> getproperty(elements[type], :vdw_radius),
    output_dots::Bool=false,
    unitcell::Union{AbstractVector,AbstractMatrix,Nothing}=nothing,
    parallel=true,
    N_SIMD::Val{N}=Val(16), # Size of SIMD blocks. Can be tuned for maximum performance.
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

                By default, sasa_particles uses the element of the atom to define the atom type,
                deducing the element from the atom name. Please verify the provided atom names.

            """))
        end
        dot_cache[type] = generate_dots(atom_radius_from_type(type), probe_radius, n_dots)
    end

    system = ParticleSystem(
        xpositions=positions(atoms),
        unitcell=unitcell,
        cutoff=2 * (maximum(atom_radius_from_type(type) for type in atom_types) + probe_radius),
        output=AtomDotMatrix(ones(Bool, n_dots, length(atoms))),
        output_name=:surface_dots,
        parallel=parallel,
    )

    # Compute sasa (function barrier for specialization on the output type)
    return _compute_sasa_particles(
        RadiiType,
        system,
        atoms,
        n_dots,
        dot_cache,
        output_dots,
        atom_type,
        atom_radius_from_type,
        probe_radius,
        N_SIMD,
    )
end

function _compute_sasa_particles(
    RadiiType::Type{<:AtomicRadiiType},
    system,
    atoms,
    n_dots,
    dot_cache,
    output_dots,
    atom_type,
    atom_radius_from_type,
    probe_radius,
    N_SIMD,
)

    pairwise!(
        (pair, surface_dots) -> update_pair_dot_exposure!(
                pair, surface_dots;
                atoms, dot_cache, atom_type, atom_radius_from_type, probe_radius, N_SIMD,
            ),
        system,
    )

    s = zeros(Float32, length(atoms))
    dots = Vector{SVector{3,Float32}}[]
    for i in eachindex(atoms)
        type_i = atom_type(atoms[i])
        s[i] = 4π * (atom_radius_from_type(atom_type(atoms[i])) + probe_radius)^2 * sum(system.surface_dots[i].exposed) / n_dots
        if output_dots
            dot_cache_i = dot_cache[type_i]
            exposed_i = system.surface_dots[i].exposed
            dots_exposed_i = SVector{3,Float32}[]
            for idot in 1:n_dots
                if exposed_i[idot]
                    x = SVector(atoms[i].x, atoms[i].y, atoms[i].z)
                    v = SVector(dot_cache_i.x[idot], dot_cache_i.y[idot], dot_cache_i.z[idot])
                    push!(dots_exposed_i, x + v)
                end
            end
            push!(dots, dots_exposed_i)
        end
    end
    return SASA{RadiiType,3,typeof(atoms)}(atoms, s, dots)
end

"""
    sasa_particles(StandardAtomicRadii, atoms; kargs...)

Alias for computing `sasa_particles` with the default standard atomistic radii (all-atom
representation of the structure, assuming the presence of hydrogen atoms). This is 
similar to calling `sasa_particles(atoms; kargs..)` omittting the first argument.

"""
function sasa_particles(::Type{StandardAtomicRadii}, atoms; kargs...)
    return _sasa_particles(
        StandardAtomicRadii, atoms; 
        atom_type=element,
        atom_radius_from_type=type -> getproperty(elements[type], :vdw_radius),
        kargs...
    )
end

"""
    sasa(s::SASA)
    sasa(s::SASA{<:AbstractVector{<:PDBTools.Atom}})
    sasa(atoms::SASA{<:AbstractVector{PDBTools.Atom}}, selection::Union{Function,String})

Given the output of `sasa_particles`, sums up contributions of atoms to compute the SASA
of the full structure, an atom, or a subset of atoms. The function can be called with only a 
`SASA` object (in which case the full SASA is returned),
or with the object and a selection, given by a function or selection string. 

# Example

```julia-repl
julia> using PDBTools

julia> prot = select(read_pdb(PDBTools.TESTPDB), "protein");

julia> at_sasa = sasa_particles(prot);

julia> sasa(at_sasa) # total sasa of prot
5389.0146f0

julia> sasa(at_sasa, "backbone") # selection string
988.7648f0

julia> sasa(at_sasa, at -> name(at) == "CA") # selection function
44.078426f0

julia> sasa(at_sasa[1]) # single atom SASA
5.467941f0
```

"""
function sasa end
sasa(p::SASA) = sum(i -> p[i], eachindex(p); init=0.f0)
sasa(p::SASA{T,N,<:AbstractVector{<:PDBTools.Atom}}, sel::Function) where {T,N} = sum(p[i] for i in eachindex(p) if sel(p.particles[i]); init=0.f0)
sasa(p::SASA{T,N,<:AbstractVector{<:PDBTools.Atom}}, sel::AbstractString) where {T,N} = sasa(p, Select(sel))

#
# Save and load SASA objects
#
function _sasa_radii_type_name(::Type{R}) where {R<:AtomicRadiiType}
    return replace(string(R), "PDBTools." => "")
end

function _sasa_radii_type(s::String)
    s == "StandardAtomicRadii" && return StandardAtomicRadii
    s == "CustomAtomicRadii" && return CustomAtomicRadii
    if isdefined(PDBTools, :CreamerUnitedAtomRadii) && s == "CreamerUnitedAtomRadii"
        return PDBTools.CreamerUnitedAtomRadii
    end
    throw(ArgumentError("Invalid SASA radii type: $s"))
end

function _atom_to_json(at::Atom)
    if at.custom !== nothing
        throw(ArgumentError("Saving SASA with atom custom fields is not supported."))
    end
    return Dict(
        "index" => at.index,
        "index_pdb" => at.index_pdb,
        "name" => String(at.name),
        "resname" => String(at.resname),
        "chain" => String(at.chain),
        "resnum" => at.resnum,
        "residue" => at.residue,
        "x" => at.x,
        "y" => at.y,
        "z" => at.z,
        "beta" => at.beta,
        "occup" => at.occup,
        "model" => at.model,
        "segname" => String(at.segname),
        "pdb_element" => String(at.pdb_element),
        "charge" => at.charge,
        "flag" => at.flag,
    )
end

function _atom_from_json(d::JSON.Object)
    return Atom{Nothing}(;
        index=Int32(d["index"]),
        index_pdb=Int32(d["index_pdb"]),
        name=String(d["name"]),
        resname=String(d["resname"]),
        chain=String(d["chain"]),
        resnum=Int32(d["resnum"]),
        residue=Int32(d["residue"]),
        x=Float32(d["x"]),
        y=Float32(d["y"]),
        z=Float32(d["z"]),
        beta=Float32(d["beta"]),
        occup=Float32(d["occup"]),
        model=Int32(d["model"]),
        segname=String(d["segname"]),
        pdb_element=String(d["pdb_element"]),
        charge=Float32(d["charge"]),
        flag=Int8(d["flag"]),
    )
end

"""
    save(filename::AbstractString, s::SASA)

Save `SASA` object data to `filename` (json format).
Load with `load(SASA, filename)`.

"""
function save(filename::AbstractString, s::SASA{R,N,<:AbstractVector{<:Atom}}) where {R<:AtomicRadiiType,N}
    dots = [
        [[Float32(x[1]), Float32(x[2]), Float32(x[3])] for x in dotset]
        for dotset in s.dots
    ]
    data = Dict(
        "radii_type" => _sasa_radii_type_name(R),
        "ndim" => N,
        "particles" => _atom_to_json.(s.particles),
        "sasa" => s.sasa,
        "dots" => dots,
    )
    open(expanduser(filename), "w") do io
        JSON.print(io, data)
    end
    return filename
end

"""
    load(SASA, filename::String)

Creates a `SASA` object from the data saved to `filename`, with the
`save(filename, s)` function.

"""
function load(::Type{SASA}, filename::AbstractString)
    data = open(expanduser(filename), "r") do io
        JSON.parse(io)
    end
    required_keys = ("radii_type", "ndim", "particles", "sasa", "dots")
    for k in required_keys
        haskey(data, k) || throw(ArgumentError("Invalid SASA file: missing key \"$k\"."))
    end
    R = _sasa_radii_type(String(data["radii_type"]))
    N = Int(data["ndim"])
    particles = _atom_from_json.(data["particles"])
    sasa_values = Float32.(data["sasa"])
    if N != 3
        throw(ArgumentError("Invalid SASA file: only 3D dots are supported, got N = $N."))
    end
    dots = Vector{Vector{SVector{3,Float32}}}(undef, length(data["dots"]))
    for i in eachindex(dots)
        raw_dotset = data["dots"][i]
        dots[i] = [SVector{3,Float32}(Float32(v[1]), Float32(v[2]), Float32(v[3])) for v in raw_dotset]
    end
    return SASA{R,3,typeof(particles)}(particles, sasa_values, dots)
end

@testitem "sasa" begin
    using PDBTools
    using ShowMethodTesting
    prot = read_pdb(PDBTools.TESTPDB, "protein")
    N = 10_030 # number of dot samples: not multiple of 16 to have a SIMD remaining

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
    at_sasa = sasa_particles(prot; n_dots=N, atom_radius_from_type=type -> vmd_radii[type])
    @test sasa(at_sasa) ≈ 5365.55029296875 rtol = 0.01
    # Accessibility of groups within the structure
    @test sasa(at_sasa, "backbone") ≈ 1130.37646484375 rtol = 0.05
    @test sasa(at_sasa, "resname GLU LYS") ≈ 797.8261108398438 rtol = 0.05
    @test sasa(at_sasa, "residue 1") ≈ 124.57905578613281 rtol = 0.05
    @test sasa(at_sasa, "residue 104") ≈ 122.50507354736328 rtol = 0.05

    # Test periodic boundary conditions
    at_sasa_no_pbc = sasa_particles(read_pdb(PDBTools.TESTNOPBC, "protein"); n_dots=N)
    uc = [107.845, 107.845, 107.845]
    at_sasa_pbc = sasa_particles(read_pdb(PDBTools.TESTPBC, "protein"); unitcell=uc, n_dots=N)
    @test sasa(at_sasa_pbc) ≈ sasa(at_sasa_no_pbc)

    # Compare with Gromacs - 2023.3 output
    # gmx sasa -s prot.pdb -o sasa_output.xvg -ndots 100000
    @test sasa(sasa_particles(prot; n_dots=N)) ≈ 100 * 53.754 rtol = 0.01
    # Isolated groups
    @test sasa(sasa_particles(select(prot, "backbone and not name O"); n_dots=N)) ≈ 100 * 55.229 rtol = 0.05
    @test sasa(sasa_particles(select(prot, "name CA"); n_dots=N)) ≈ 100 * 58.630 rtol = 0.05
    @test sasa(sasa_particles(select(prot, "sidechain and not element H"); n_dots=N)) ≈ 100 * 69.024 rtol = 0.05

    # Test non-contiguous indexing with general selections
    at_sasa = sasa_particles(select(prot, "name CA"); n_dots=N)
    @test sasa(at_sasa) ≈ 5863.24f0
    @test sasa(at_sasa, "resname THR") ≈ 322.17105f0
    @test sasa(at_sasa, at -> resname(at) == "THR") ≈ 322.17105f0

    # Test parallelization
    @test sasa(sasa_particles(prot; parallel=false)) ≈ sasa(sasa_particles(prot; parallel=true))

    # Test the output of dots
    at_sasa = sasa_particles(select(prot, "name CA"); n_dots=20, output_dots=true)
    @test sasa(at_sasa) ≈ 5856.9966f0
    @test sum(length(v) for v in at_sasa.dots) == 970

    # Test save/load for SASA with StandardAtomicRadii
    tmp = tempname() * ".json"
    save(tmp, at_sasa)
    at_sasa_load = load(SASA, tmp)
    @test at_sasa_load isa PDBTools.SASA{StandardAtomicRadii,3,Vector{Atom{Nothing}}}
    @test at_sasa_load.particles == at_sasa.particles
    @test at_sasa_load.sasa == at_sasa.sasa
    @test at_sasa_load.dots == at_sasa.dots
    rm(tmp; force=true)

    # Test save/load for SASA with CreamerUnitedAtomRadii
    at_sasa_creamer = sasa_particles(PDBTools.CreamerUnitedAtomRadii, prot; n_dots=20, output_dots=true)
    save(tmp, at_sasa_creamer)
    at_sasa_creamer_load = load(PDBTools.SASA, tmp)
    @test at_sasa_creamer_load isa PDBTools.SASA{PDBTools.CreamerUnitedAtomRadii,3,Vector{Atom{Nothing}}}
    @test at_sasa_creamer_load.particles == at_sasa_creamer.particles
    @test at_sasa_creamer_load.sasa == at_sasa_creamer.sasa
    @test at_sasa_creamer_load.dots == at_sasa_creamer.dots 
    rm(tmp; force=true)

    # Unsupported save input: custom atom fields are not serializable
    atoms_custom = [Atom(; custom=Dict(:x => 1), name="C", pdb_element="C")]
    at_sasa_custom = sasa_particles(atoms_custom; atom_type=at -> "X", atom_radius_from_type=type -> 1.5f0, n_dots=20)
    @test_throws "custom fields is not supported" save(tmp, at_sasa_custom)

    # Unsupported load inputs: malformed or invalid serialized files
    open(tmp, "w") do io
        write(io, "{\"radii_type\":\"StandardAtomicRadii\"}")
    end
    @test_throws "missing key \"ndim\"" load(PDBTools.SASA, tmp)

    open(tmp, "w") do io
        write(io, "{\"radii_type\":\"InvalidRadii\",\"ndim\":3,\"particles\":[],\"sasa\":[],\"dots\":[]}")
    end
    @test_throws "Invalid SASA radii type" load(PDBTools.SASA, tmp)

    open(tmp, "w") do io
        write(io, "{\"radii_type\":\"StandardAtomicRadii\",\"ndim\":2,\"particles\":[],\"sasa\":[],\"dots\":[]}")
    end
    @test_throws "only 3D dots are supported" load(PDBTools.SASA, tmp)
    rm(tmp; force=true)

    # Test show method
    @test parse_show(at_sasa; repl=Dict(r"PDBTools." => "")) ≈ """
            SASA{StandardAtomicRadii, 3, Vector{Atom{Nothing}}}
                Number of particles: 104
                Total SASA: 5856.9966
                Output of dots: true 
        """

    # Test errors
    @test_throws "probe_radius too small" sasa_particles(prot; probe_radius=-2.0)
    @test_throws "n_dots incorrectly" sasa_particles(prot; n_dots=0)
    prot[1].name = "Ti"
    prot[1].pdb_element = "Ti"
    @test_throws "Ti does not have" sasa_particles(prot)

end

