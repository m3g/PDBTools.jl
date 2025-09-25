import MolSimToolkitShared: dihedral

"""
    dihedral(at1::Atom, at2::Atom, at3::Atom, at4::Atom)

Computes the dihderal angle given four atoms of type PDBTools.Atom.

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
    dihedral(coor(at1), coor(at2), coor(at3), coor(at4))

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
    print(chomp("""
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
    pdb = read_pdb(PDBTools.TESTPDB)
    @test_throws ArgumentError Ramachandran(pdb)
    prot = filter(isprotein, pdb)
    pp = copy.(prot)
    popat!(pp, 15) # delete CA of 2nd residue
    @test_throws ArgumentError Ramachandran(pp)
    ram = Ramachandran(prot)
    @test sum(ram.phi) ≈ -7024.354f0
    @test sum(ram.psi) ≈ 3501.747f0
end
