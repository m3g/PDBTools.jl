import MolSimToolkitShared: dihedral

"""
    dihedral(at1::Atom, at2::Atom, at3::Atom, at4::Atom)

Computes the dihedral angle given four atoms of type PDBTools.Atom.

# Example

```jldoctest
julia> using PDBTools

julia> pdb = read_pdb(PDBTools.TESTPDB);

julia> C1 = pdb[11]; N2 = pdb[13]; CA2 = pdb[15]; C2 = pdb[22];

julia> phi = dihedral(C1, N2, CA2, C2) 
-36.70359f0
```

"""
dihedral(at1::Atom, at2::Atom, at3::Atom, at4::Atom) =
    dihedral(position(at1), position(at2), position(at3), position(at4))

@testitem "dihedral" begin
    using PDBTools
    pdb = read_pdb(PDBTools.TESTPDB);
    N1 = pdb[1]; CA1 = pdb[5]; C1 = pdb[11]; N2 = pdb[13];
    @test dihedral(N1, CA1, C1, N2) ≈ 64.07296f0
    C1 = pdb[11]; N2 = pdb[13]; CA2 = pdb[15]; C2 = pdb[22];
    @test dihedral(C1, N2, CA2, C2) ≈ -36.70359f0
end

"""
    Ramachandran(prot::AbstractVector{<:PDBTools.Atom})
    Ramachandran # type

The `Ramachandran` function receives a vector of atoms of a protein and
and returns a `Ramachandran` object, with two fields `phi` and `psi`, containing
the lists of corresponding angles, that is:

- phi: C(-1) - N - CA - C
- psi: N - CA - C - N(+1)

If any of the above atoms is missing, the function errors. The residues are expected
to belong to a single chain and consecutive. 

The resulting `Ramachandran` object can be plotted with the `Plots.scatter` function.

# Example

```jldoctest
julia> using PDBTools

julia> prot = read_pdb(PDBTools.TESTPDB, "protein");

julia> ram = Ramachandran(prot)
Ramachandran data: phi, psi vectors with 102 angles.

```

"""
@kwdef struct Ramachandran
    phi::Vector{Float32} = Float32[]
    psi::Vector{Float32} = Float32[]
end

function Base.show(io::IO, ram::Ramachandran)
    print(io, chomp("""
    Ramachandran data: phi, psi vectors with $(length(ram.phi)) angles.
    """))
end

function _fetch_atom(atname, residue)
    iat = findfirst(at -> name(at) == atname, residue)
    if isnothing(iat)
        throw(ArgumentError("""\n
            Could not find atom $atname in residue $(resname(residue))$(resnum(residue)))
        """))
    end
    return residue[iat]
end

function Ramachandran(prot::AbstractVector{<:Atom})
    ram = Ramachandran()
    residues = collect(eachresidue(prot))
    for ires in firstindex(residues)+1:lastindex(residues)-1
        rlast = residues[ires-1]
        res = residues[ires]
        rnext = residues[ires+1]
        for r in (rlast, res, rnext)
            if !isprotein(r)
                throw(ArgumentError("""\n
                    Residue $(resname(r)) in the provided array of atoms is not a protein residue.
                    Ramachandran plot cannot be computed.
                """))
            end
        end
        C0 = _fetch_atom("C", rlast)
        N = _fetch_atom("N", res)
        CA = _fetch_atom("CA", res)
        C = _fetch_atom("C", res)
        push!(ram.phi, dihedral(C0, N, CA, C))
        N2 = _fetch_atom("N", rnext)
        N2 = rnext[findfirst(at -> name(at) == "N", rnext)]
        push!(ram.psi, dihedral(N, CA, C, N2))
    end
    return ram
end

@testitem "Ramachandran" begin
    using PDBTools
    using ShowMethodTesting
    pdb = read_pdb(PDBTools.TESTPDB)
    @test_throws ArgumentError Ramachandran(pdb)
    prot = filter(isprotein, pdb)
    pp = copy.(prot)
    popat!(pp, 15) # delete CA of 2nd residue
    @test_throws ArgumentError Ramachandran(pp)
    ram = Ramachandran(prot)
    @test sum(ram.phi) ≈ -7024.354f0
    @test sum(ram.psi) ≈ 3501.747f0
    @test parse_show(ram) ≈ "Ramachandran data: phi, psi vectors with 102 angles."
end

function _validate_protein_chain(atoms::AbstractVector{<:Atom})
    isempty(atoms) && throw(ArgumentError("The vector of atoms is empty."))
    if !Base.all(isprotein, atoms)
        throw(ArgumentError("""\n
            All atoms must be protein atoms to set backbone dihedral angles.
        """))
    end
    chains = unique(chain.(atoms))
    if length(chains) > 1
        throw(ArgumentError("""\n
            The atoms must belong to a single chain to set backbone dihedral angles. Got chains: $(chains).
        """))
    end
    return collect(eachresidue(atoms))
end

function _angle_in_degrees(angle::Real, unit::AbstractString)
    if unit == "rad"
        return rad2deg(angle)
    elseif unit == "deg"
        return float(angle)
    else
        throw(ArgumentError("""\n
            Invalid unit: "$unit". Use "rad" or "deg".
        """))
    end
end

# Rigidly rotates `moving_atoms` by `angle_deg` around the line defined by `origin` and
# `axis` (which need not be normalized nor coincide with either dihedral pivot atom, since
# a rotation about a line does not depend on the reference point chosen along it).
function _rotate_downstream!(moving_atoms, origin::AbstractVector, axis::AbstractVector, angle_deg::Real)
    k = axis / norm(axis)
    sinθ, cosθ = sincosd(angle_deg)
    for at in moving_atoms
        v = position(at) - origin
        v_rot = v * cosθ + cross(k, v) * sinθ + k * (dot(k, v) * (1 - cosθ))
        set_position!(at, origin + v_rot)
    end
    return nothing
end

"""
    set_phi!(atoms::AbstractVector{<:Atom}, phi_vec::AbstractVector{<:Real}; unit="rad")

Sets the phi backbone dihedral angles (C(-1) - N - CA - C) of a protein, given a vector of
target angles. `phi_vec` must have `length(residues) - 1` elements, one for every residue
but the first (in sequential order), where `residues = collect(eachresidue(atoms))`.

Angles are interpreted as radians by default (`unit="rad"`); pass `unit="deg"` to provide
degrees instead.

Setting the phi angle of a residue rigidly rotates every atom attached to its `CA`, other
than `N` - that is, its side chain, its `C` and `O` (and `OXT`, if terminal) atoms, and all
atoms of the following residues - around the `N`-`CA` bond, so that the new phi angle
matches the requested value. Rotations are applied residue by residue, from the first to
the last, so that previously-set phi angles are not disturbed by later ones.

`atoms` must belong to a single protein chain, and contain the `N`, `CA`, and `C` atoms
needed to define every phi angle to be set, or the function throws an `ArgumentError`.

# Example

```jldoctest
julia> using PDBTools

julia> prot = read_pdb(PDBTools.TESTPDB, "protein");

julia> nres = length(collect(eachresidue(prot)));

julia> set_phi!(prot, fill(-60.0, nres - 1); unit="deg");

julia> ram = Ramachandran(prot);

julia> all(≈(-60.0; atol=1e-3), ram.phi)
true

```

"""
function set_phi!(atoms::AbstractVector{<:Atom}, phi_vec::AbstractVector{<:Real}; unit="rad")
    residues = _validate_protein_chain(atoms)
    n = length(residues)
    n >= 2 || throw(ArgumentError("At least two residues are required to define a phi angle."))
    if length(phi_vec) != n - 1
        throw(ArgumentError("""\n
            Expected $(n - 1) phi angles (one for every residue but the first), got $(length(phi_vec)).
        """))
    end
    for i in 2:n
        rprev = residues[i-1]
        res = residues[i]
        C0 = _fetch_atom("C", rprev)
        N = _fetch_atom("N", res)
        CA = _fetch_atom("CA", res)
        C = _fetch_atom("C", res)
        current = dihedral(C0, N, CA, C)
        target = _angle_in_degrees(phi_vec[i-1], unit)
        moving_atoms = Iterators.flatten((
            (at for at in res if at !== N),
            (at for r in @view(residues[i+1:end]) for at in r),
        ))
        _rotate_downstream!(moving_atoms, position(N), position(CA) - position(N), target - current)
    end
    return atoms
end

"""
    set_psi!(atoms::AbstractVector{<:Atom}, psi_vec::AbstractVector{<:Real}; unit="rad")

Sets the psi backbone dihedral angles (N - CA - C - N(+1)) of a protein, given a vector of
target angles. `psi_vec` must have `length(residues) - 1` elements, one for every residue
but the last (in sequential order), where `residues = collect(eachresidue(atoms))`.

Angles are interpreted as radians by default (`unit="rad"`); pass `unit="deg"` to provide
degrees instead.

Setting the psi angle of a residue rigidly rotates every atom attached to its `C`, other
than `CA` - that is, its `O` atom and all atoms of the following residues - around the
`CA`-`C` bond, so that the new psi angle matches the requested value. Rotations are applied
residue by residue, from the first to the last, so that previously-set psi (and phi) angles
are not disturbed by later ones.

`atoms` must belong to a single protein chain, and contain the `N`, `CA`, `C`, and `O`
atoms needed to define every psi angle to be set, or the function throws an `ArgumentError`.

# Example

```jldoctest
julia> using PDBTools

julia> prot = read_pdb(PDBTools.TESTPDB, "protein");

julia> nres = length(collect(eachresidue(prot)));

julia> set_psi!(prot, fill(140.0, nres - 1); unit="deg");

julia> ram = Ramachandran(prot);

julia> all(≈(140.0; atol=1e-3), ram.psi)
true

```

"""
function set_psi!(atoms::AbstractVector{<:Atom}, psi_vec::AbstractVector{<:Real}; unit="rad")
    residues = _validate_protein_chain(atoms)
    n = length(residues)
    n >= 2 || throw(ArgumentError("At least two residues are required to define a psi angle."))
    if length(psi_vec) != n - 1
        throw(ArgumentError("""\n
            Expected $(n - 1) psi angles (one for every residue but the last), got $(length(psi_vec)).
        """))
    end
    for i in 1:n-1
        res = residues[i]
        rnext = residues[i+1]
        N = _fetch_atom("N", res)
        CA = _fetch_atom("CA", res)
        C = _fetch_atom("C", res)
        O = _fetch_atom("O", res)
        N2 = _fetch_atom("N", rnext)
        current = dihedral(N, CA, C, N2)
        target = _angle_in_degrees(psi_vec[i], unit)
        moving_atoms = Iterators.flatten((
            (O,),
            (at for r in @view(residues[i+1:end]) for at in r),
        ))
        _rotate_downstream!(moving_atoms, position(CA), position(C) - position(CA), target - current)
    end
    return atoms
end



"""
    extended_chain!(atoms::AbstractVector{<:Atom})

Sets every phi and psi backbone dihedral angle of `atoms` to 180°, producing the fully
extended (all-trans, "C5") conformation, in place.

`atoms` must satisfy the same requirements as `set_phi!`/`set_psi!` (single protein chain,
required backbone atoms present).

"""
function extended_chain!(atoms::AbstractVector{<:Atom})
    nres = length(collect(eachresidue(atoms)))
    theta = fill(180.0, nres - 1)
    set_phi!(atoms, theta; unit="deg")
    set_psi!(atoms, theta; unit="deg")
    return atoms
end

"""
    extended_chain(atoms::AbstractVector{<:Atom})

Non-mutating version of `extended_chain!`: returns a new vector of atoms, a copy of `atoms`
with every phi and psi backbone dihedral angle set to 180° (the fully extended conformation).

"""
function extended_chain(atoms::AbstractVector{<:Atom})
    ext = copy.(atoms)
    extended_chain!(ext)
    return ext
end

@testitem "extended_chain and extended_chain!" begin
    using PDBTools

    prot = read_pdb(PDBTools.TESTPDB, "protein")
    ext = extended_chain(prot)

    # Non-mutating: the original is untouched
    orig = read_pdb(PDBTools.TESTPDB, "protein")
    for (at1, at2) in zip(prot, orig)
        @test position(at1) ≈ position(at2)
    end

    ram = Ramachandran(ext)
    @test all(x -> isapprox(x, 180.0; atol=1e-2) || isapprox(abs(x), 180.0; atol=1e-2), ram.phi)
    @test all(x -> isapprox(x, 180.0; atol=1e-2) || isapprox(abs(x), 180.0; atol=1e-2), ram.psi)

    # Mutating version changes the input in place
    prot2 = read_pdb(PDBTools.TESTPDB, "protein")
    extended_chain!(prot2)
    ram2 = Ramachandran(prot2)
    @test all(x -> isapprox(x, 180.0; atol=1e-2) || isapprox(abs(x), 180.0; atol=1e-2), ram2.phi)
    @test all(x -> isapprox(x, 180.0; atol=1e-2) || isapprox(abs(x), 180.0; atol=1e-2), ram2.psi)
end

@testitem "set_phi! and set_psi!" begin
    using PDBTools

    prot = read_pdb(PDBTools.TESTPDB, "protein")
    nres = length(collect(eachresidue(prot)))

    # Setting all phi angles to a fixed value reproduces that value everywhere
    set_phi!(prot, fill(-60.0, nres - 1); unit="deg")
    ram = Ramachandran(prot)
    @test all(x -> isapprox(x, -60.0; atol=1e-2), ram.phi)

    # Setting all psi angles to a fixed value reproduces that value everywhere,
    # and does not disturb the phi angles set above
    set_psi!(prot, fill(140.0, nres - 1); unit="deg")
    ram = Ramachandran(prot)
    @test all(x -> isapprox(x, -60.0; atol=1e-2), ram.phi)
    @test all(x -> isapprox(x, 140.0; atol=1e-2), ram.psi)

    # Round-trip: extracting a structure's own (arbitrary) dihedrals and re-applying
    # them onto an identical, independent copy must reproduce its coordinates. Since
    # `prot` at this point is the non-native (-60, 140) conformation built above, this
    # is a real check of the rotated-fragment logic, not a trivial zero-angle no-op.
    residues = collect(eachresidue(prot))
    phi_native = Float64[]
    psi_native = Float64[]
    for i in 2:nres
        rprev, res = residues[i-1], residues[i]
        push!(phi_native, dihedral(
            PDBTools._fetch_atom("C", rprev), PDBTools._fetch_atom("N", res),
            PDBTools._fetch_atom("CA", res), PDBTools._fetch_atom("C", res),
        ))
    end
    for i in 1:nres-1
        res, rnext = residues[i], residues[i+1]
        push!(psi_native, dihedral(
            PDBTools._fetch_atom("N", res), PDBTools._fetch_atom("CA", res),
            PDBTools._fetch_atom("C", res), PDBTools._fetch_atom("N", rnext),
        ))
    end
    prot_copy = read_pdb(PDBTools.TESTPDB, "protein")
    set_phi!(prot_copy, fill(-60.0, nres - 1); unit="deg")
    set_psi!(prot_copy, fill(140.0, nres - 1); unit="deg")
    set_phi!(prot_copy, phi_native; unit="deg")
    set_psi!(prot_copy, psi_native; unit="deg")
    for (at1, at2) in zip(prot_copy, prot)
        @test position(at1) ≈ position(at2) atol = 1e-1
    end

    # radians is the default unit
    prot3 = read_pdb(PDBTools.TESTPDB, "protein")
    set_phi!(prot3, fill(-1.0, nres - 1))
    ram3 = Ramachandran(prot3)
    @test all(x -> isapprox(x, rad2deg(-1.0); atol=1e-2), ram3.phi)

    # Errors
    @test_throws ArgumentError set_phi!(prot, fill(0.0, nres)) # wrong length
    @test_throws ArgumentError set_psi!(prot, fill(0.0, nres)) # wrong length
    @test_throws ArgumentError set_phi!(prot, fill(0.0, nres - 1); unit="degrees") # invalid unit

    pdb = read_pdb(PDBTools.TESTPDB) # not a protein-only selection
    @test_throws ArgumentError set_phi!(pdb, fill(0.0, length(collect(eachresidue(pdb))) - 1))

    pp = copy.(prot)
    popat!(pp, findfirst(at -> name(at) == "CA" && resnum(at) == resnum(residues[2][1]), pp))
    @test_throws ArgumentError set_phi!(pp, fill(0.0, nres - 1))

    # multiple chains
    dimer = read_pdb(PDBTools.TESTPDB, "protein")
    for at in @view(dimer[(length(dimer) ÷ 2 + 1):end])
        at.chain = "Z"
    end
    @test_throws ArgumentError set_phi!(dimer, fill(0.0, nres - 1))
end