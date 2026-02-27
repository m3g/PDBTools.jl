
using ProteinSecondaryStructures: 
    ProteinSecondaryStructures,
    stride_run,
    dssp_run,
    ss_composition,
    ss_name,
    ss_number,
    SSData

export stride_run, dssp_run, ss_composition, ss_name, ss_number, SSData

#= 

This function replaces the residue names with their standard three-letter codes,
and checks for the presence of backbone atoms (N, CA, C, O) in each residue.

The absence of any backbone atom will trigger a warning message, because with 
missing backbone atoms, secondary structure assignment tools like STRIDE and DSSP
may not function correctly.

=# 
function _set_pdb(ats::AbstractVector{<:PDBTools.Atom})
    ats_new = copy.(select(ats, at -> isprotein(at) && element(at) != "H"))
    for r in eachresidue(ats_new)
        N = false
        CA = false
        C = false
        O = false
        for at in r
            at.resname = threeletter(resname(at))
            name(at) == "N" && (N = true)
            name(at) == "CA" && (CA = true)
            name(at) == "C" && (C = true)
            name(at) == "O" && (O = true)
        end
        N || @warn("Residue $(resname(r))$(resnum(r)) - chain $(chain(r)) is missing N backbone atom.")
        CA || @warn("Residue $(resname(r))$(resnum(r)) - chain $(chain(r)) is missing CA backbone atom.")
        C || @warn("Residue $(resname(r))$(resnum(r)) - chain $(chain(r)) is missing C backbone atom.")
        O || @warn("Residue $(resname(r))$(resnum(r)) - chain $(chain(r)) is missing O backbone atom.")
    end
    return ats_new
end

"""
    stride_run(atoms::AbstractVector{<:PDBTools.Atom})

Run STRIDE secondary structure assignment on the provided array of atoms.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.TESTPDB, "protein");

julia> atoms[1458].name = "O"; # Terminal residue has non-standard atom name

julia> ss = stride_run(atoms)
104-element Vector{SSData}:
 SSData("ALA", "A", 1, "C", 360.0, 64.07)
 SSData("CYS", "A", 2, "T", -36.7, 125.51)
 SSData("ASP", "A", 3, "T", -125.56, -10.38)
 ⋮
 SSData("CYS", "A", 103, "C", -56.48, -168.79)
 SSData("THR", "A", 104, "C", -110.75, 360.0)
```

"""
function ProteinSecondaryStructures.stride_run(atoms::AbstractVector{<:PDBTools.Atom})
    ats_new = _set_pdb(atoms)
    # Possibly remove non-standard chain names coming from CIF data
    non_standard_chain = findfirst(c -> length(c) > 1, chain(at) for at in ats_new)
    if !isnothing(non_standard_chain)
        original_chains = chain.(eachresidue(ats_new))
        c = String15("A")
        iat = 0
        ires = 1
        for r in eachresidue(ats_new)
            if ires == 1
                iat += length(r)
                ires += 1
                continue
            end
            if original_chains[ires] != original_chains[ires-1]
                c = c == String15("A") ? String15("B") : String15("A")
            end
            for at in r
                iat += 1
                ats_new[iat].chain = c
            end
            ires += 1
        end
    end
    tmppdb = tempname() * ".pdb"
    write_pdb(tmppdb, ats_new)
    ss = stride_run(tmppdb; adjust_pdb=true)
    if !isnothing(non_standard_chain)
        if length(ss) != length(original_chains)
            @warn("""\n
                Mapping of stride output residues failed relative to original residue names,
                which do not follow the PDB single-letter standard.

                    Number of original residues: $(length(original_chains))
                    Number of resulting secondary structure data elements: $(length(ss))

                The chain identifiers of the secondary structure elements might be incorrect.

            """)
        end
        is = 0
        for s in ss
            is += 1
            ss[is] = SSData(s.resname, original_chains[is], s.resnum, s.sscode, s.phi, s.psi)
        end
    end
    rm(tmppdb; force=true)
    return ss
end

"""
    dssp_run(atoms::AbstractVector{<:PDBTools.Atom})

Run DSSP secondary structure assignment on the provided array of atoms.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.TESTPDB, "protein");

julia> atoms[1458].name = "O"; # Terminal residue has non-standard atom name

julia> ss = dssp_run(atoms)
104-element Vector{SSData}:
 SSData("ALA", "A", 1, " ", 0.0, 64.1)
 SSData("CYS", "A", 2, " ", -36.7, 125.5)
 SSData("ASP", "A", 3, " ", -125.6, -10.4)
 ⋮
 SSData("CYS", "A", 103, " ", -56.5, -168.8)
 SSData("THR", "A", 104, " ", -110.7, 0.0)
```

"""
function ProteinSecondaryStructures.dssp_run(atoms::AbstractVector{<:PDBTools.Atom})
    ats_new = _set_pdb(atoms)
    tmpcif = tempname() * ".cif"
    write_mmcif(tmpcif, ats_new)
    ss = dssp_run(tmpcif)
    rm(tmpcif; force=true)
    return ss
end

@testitem "secondary structure assignment" begin
    using PDBTools
    ats = read_pdb(PDBTools.TESTPDB, "protein")
    ss_stride = stride_run(ats)
    ss_dssp = dssp_run(ats)
    @test length(ss_stride) == 104
    @test length(ss_dssp) == 103
    @test ss_composition(ss_stride) == Dict("310 helix" => 0, "bend" => 0, "turn" => 36, "beta bridge" => 4, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 21, "alpha helix" => 17, "loop" => 0, "coil" => 26)
    @test ss_composition(ss_dssp) == Dict("310 helix" => 0, "bend" => 18, "turn" => 9, "beta bridge" => 3, "kappa helix" => 2, "pi helix" => 0, "beta strand" => 15, "alpha helix" => 17, "loop" => 39, "coil" => 0) 
end