
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
function _set_pdb(ats::AbstractVector{<:Atom})
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

function _chain_map(i)
    c = if i <= 94
        c = ('\x21':'\x7e')[i]
    else
        throw(ArgumentError(("""\n
            Number of chains greater than 94 not supported by stride. Use dssp_run instead.

        """)))
    end
    return String15(string(c))
end

#=
DSSP and STRIDE skip residues they cannot classify (e.g. residues with missing backbone
atoms). Here, the skipped residues are added back to the output, with the code of the 
program for residues without a secondary structure (`sscode`, "loop" in DSSP and 
"coil" in STRIDE) and undefined dihedral angles (`undefined_angle`), such that the 
output has one element per residue of `ats`, in the same order.
=#
function _add_missing_residues(ss::AbstractVector{SSData}, ats::AbstractVector{<:Atom}; 
    sscode::String, undefined_angle::Float64, method::String,
)
    ss_full = SSData[]
    sizehint!(ss_full, length(ss))
    # The programs may rename the chains (e.g. blank chains are reported as "?" by DSSP):
    # the chain names of the output are mapped to those of the input as they are matched.
    chain_map = Dict{String,String}()
    j = 1
    for r in eachresidue(ats)
        if j <= length(ss) && ss[j].resnum == resnum(r) && ss[j].resname == resname(r) &&
           get(chain_map, ss[j].chain, chain(r)) == chain(r)
            chain_map[ss[j].chain] = chain(r)
            push!(ss_full, ss[j])
            j += 1
        else
            @warn("Residue $(resname(r))$(resnum(r)) - chain $(chain(r)) was not assigned by $method: set as $(ss_name(sscode)).")
            push!(ss_full, SSData(resname(r), chain(r), resnum(r), sscode, undefined_angle, undefined_angle))
        end
    end
    if j <= length(ss)
        throw(ArgumentError("""\n
            The output of $method could not be matched to the residues of the structure.
            Element $j of the output: $(ss[j]) was not found in the input residues.

        """))
    end
    return ss_full
end


"""
    stride_run(atoms::AbstractVector{<:PDBTools.Atom})

Run STRIDE secondary structure assignment on the provided array of atoms.

The output has one element per protein residue of `atoms`, in the order of the residues.
Residues that STRIDE cannot classify (e.g. residues with missing backbone atoms) are
assigned as coil (`"C"`), with undefined dihedral angles (`360.0`), and a warning is shown.

# Example

```jldoctest
julia> using PDBTools, LazyArtifacts

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

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
function ProteinSecondaryStructures.stride_run(atoms::AbstractVector{<:Atom})
    ats_new = _set_pdb(atoms)
    # Possibly remove non-standard chain names coming from CIF data
    non_standard_chain = findfirst(c -> length(c) > 1, chain(at) for at in ats_new)
    if !isnothing(non_standard_chain)
        original_chains = chain.(eachresidue(ats_new))
        i = 0
        for chain in eachchain(ats_new)
            i += 1
            for at in chain
                at.chain = _chain_map(i)
            end
        end
    end
    tmppdb = tempname() * ".pdb"
    write_pdb(tmppdb, ats_new)
    ss = stride_run(tmppdb; adjust_pdb=true)
    rm(tmppdb; force=true)
    ss = _add_missing_residues(ss, ats_new; sscode="C", undefined_angle=360.0, method="STRIDE")
    # restore the original chain names: ss has one element per residue of ats_new
    if !isnothing(non_standard_chain)
        for (is, s) in enumerate(ss)
            ss[is] = SSData(s.resname, original_chains[is], s.resnum, s.sscode, s.phi, s.psi)
        end
    end
    return ss
end

"""
    dssp_run(atoms::AbstractVector{<:PDBTools.Atom})

Run DSSP secondary structure assignment on the provided array of atoms.

The output has one element per protein residue of `atoms`, in the order of the residues.
Residues that DSSP cannot classify (e.g. residues with missing backbone atoms) are
assigned as loop (`" "`), with undefined dihedral angles (`0.0`), and a warning is shown.
DSSP considers only the first model of the structure.

# Example

```jldoctest
julia> using PDBTools, LazyArtifacts

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

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
function ProteinSecondaryStructures.dssp_run(atoms::AbstractVector{<:Atom})
    ats_new = _set_pdb(atoms)
    # DSSP only considers the first model of the structure
    first_model = isempty(ats_new) ? 0 : model(first(ats_new))
    if any(at -> model(at) != first_model, ats_new)
        @warn("The structure has more than one model: DSSP assigns the secondary structure of the first model only.")
        ats_new = filter(at -> model(at) == first_model, ats_new)
    end
    tmpcif = tempname() * ".cif"
    write_mmcif(tmpcif, ats_new)
    ss = dssp_run(tmpcif)
    rm(tmpcif; force=true)
    return _add_missing_residues(ss, ats_new; sscode=" ", undefined_angle=0.0, method="DSSP")
end

@testitem "secondary structure assignment" begin
    using LazyArtifacts
    using PDBTools
    ats = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    ss_stride = stride_run(ats)
    ss_dssp = dssp_run(ats)
    @test length(ss_stride) == 104
    @test length(ss_dssp) == 104
    @test ss_composition(ss_stride) == Dict("310 helix" => 0, "bend" => 0, "turn" => 36, "beta bridge" => 4, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 21, "alpha helix" => 17, "loop" => 0, "coil" => 26)
    @test ss_composition(ss_dssp) == Dict("310 helix" => 0, "bend" => 18, "turn" => 9, "beta bridge" => 3, "kappa helix" => 2, "pi helix" => 0, "beta strand" => 15, "alpha helix" => 17, "loop" => 40, "coil" => 0)
    # THR104 has no O atom (it is named OT1): DSSP skips it, and it is added back as loop
    @test_logs (:warn, r"THR104 - chain A was not assigned by DSSP") match_mode = :any dssp_run(ats)
    @test ss_dssp[end] == SSData("THR", "A", 104, " ", 0.0, 0.0)
    # blank chain names (reported as "?" by DSSP)
    ats_blank = deepcopy(ats)
    foreach(at -> at.chain = " ", ats_blank)
    @test length(dssp_run(ats_blank)) == 104
    @test length(stride_run(ats_blank)) == 104
    # STRIDE also skips a residue without backbone atoms: added back as coil
    ats_noCA = filter(at -> !(resnum(at) == 50 && name(at) in ("N", "CA", "C", "O")), ats)
    ss = stride_run(ats_noCA)
    @test length(ss) == 104
    @test ss[50].resnum == 50
    ss = dssp_run(ats_noCA)
    @test length(ss) == 104
    @test ss[50].sscode == " "
    # multiple models: DSSP considers the first model only
    nmr = read_mmcif(PDBTools.test_file("8S8N.cif"), "protein")
    @test_logs (:warn, r"more than one model") match_mode = :any dssp_run(nmr)
    @test length(dssp_run(nmr)) == length(eachresidue(select(nmr, at -> model(at) == 1)))
    p = read_mmcif(PDBTools.test_file("6CO8.cif"))
    ss_dssp = dssp_run(p)
    @test length(ss_dssp) == 501
    @test ss_dssp[1].chain == "A_1"
    ss_stride = stride_run(p)
    @test length(ss_stride) == 501
    @test ss_stride[1].chain == "A_1"
    @test_throws "greater than 94" PDBTools._chain_map(95)
end