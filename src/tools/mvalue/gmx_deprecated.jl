#
# Deprecated: SASA per residue type computed with GROMACS (`gmx sasa`).
# Use `delta_sasa_per_restype`, which uses the internal SASA implementation.
#

#=
    read_gmx_delta_sasa_per_restype_values(filename::AbstractString, n)

Reads the output of `gmx sasa` and returns the SASA values.
`n` is the number of surfaces calculated (1 for BB only, 2 for SC and BB, for example).

=#
function read_gmx_delta_sasa_per_restype_values(filename::AbstractString, n)
    local sasa_values
    open(filename, "r") do io
        for line in eachline(io)
            if startswith(line, r"@|#")
                continue
            end
            sasa_values = parse.(Float64, split(line)[3:2+n])
        end
    end
    return sasa_values
end

"""
    gmx_delta_sasa_per_restype(; native_pdb::AbstractString, desnat_pdb::AbstractString)

Calculates the change in solvent accessible surface area (SASA) upon denaturation for each amino acid type
using GROMACS. Returns a dictionary that can be directly used as input to the `mvalue_delta_sasa` function.

!!! warning "Deprecated"
    This function is deprecated and will be removed in a future breaking release.
    Use `delta_sasa_per_restype`, which computes the SASA with the internal implementation of PDBTools.jl
    and does not require GROMACS.

!!! note
    This function requires GROMACS (`gmx sasa` executable) to be installed and accessible from the command line.
    The path to the `gmx` executable can be provided with the `gmx` keyword.

# Arguments

- `native_pdb::AbstractString`: Path to the PDB file of the native protein structure.
- `desnat_pdb::AbstractString`: Path to the PDB file of the denatured protein structure.

# Optional arguments 

- `gmx`: the path to the `gmx` GROMACS exectuable (by default it expects `gmx` to be on the path).
- `n_dots::Int`: sets the precision of the SASA grid (greater is better).
- `backbone::Function = at -> name(at) in ("N", "CA", "C", "O")`: Define what is a backbone atom.
- `sidechain::Function = at -> !(name(at) in ("N", "CA", "C", "O"))`: Define what is a sidechain atom.
- `ignore_hydrogen::Bool=true`: By default, ignore all hydrogen atoms.

# Returns

A dictionary where each key is an amino acid three-letter code (e.g., "ALA", "PHE"), and the value
is another dictionary with two keys: `:sc` for side chain SASA values and `:bb` for backbone SASA values.
Each of these keys maps to a tuple containing a single Float64 value representing the change in SASA upon denaturation in Å².

"""
function gmx_delta_sasa_per_restype(;
    native_pdb::AbstractString,
    desnat_pdb::AbstractString,
    n_dots::Int=500,
    ignore_hydrogen::Bool=true,
    gmx="gmx",
    backbone::Function=at -> name(at) in ("N", "CA", "C", "O"),
    sidechain::Function=at -> !(name(at) in ("N", "CA", "C", "O")),
)
    Base.depwarn(
        "`gmx_delta_sasa_per_restype` is deprecated and will be removed in a future breaking release. Use `delta_sasa_per_restype` instead.",
        :gmx_delta_sasa_per_restype,
    )
    sasas = OrderedDict{String,OrderedDict{Symbol,Float64}}()
    p = read_pdb(native_pdb, "protein")
    for rname in unique(resname.(eachresidue(p)))
        sasa_bb_native, sasa_sc_native = 100 .* gmx_sasa_per_restype(native_pdb, rname; gmx, n_dots, backbone, sidechain, ignore_hydrogen) # returns in nm^2
        sasa_bb_desnat, sasa_sc_desnat = 100 .* gmx_sasa_per_restype(desnat_pdb, rname; gmx, n_dots, backbone, sidechain, ignore_hydrogen)
        sasas[rname] = OrderedDict(:sc => sasa_sc_desnat - sasa_sc_native, :bb => sasa_bb_desnat - sasa_bb_native)
    end
    return sasas # returns in Å^2
end

#
# Runs gmx sasa for a single residue type in a given PDB file.
#
function gmx_sasa_per_restype(
    pdbname, rname;
    gmx="gmx",
    n_dots::Int=500,
    backbone::Function,
    sidechain::Function,
    ignore_hydrogen::Bool=true,
)
    keepH(at) = ignore_hydrogen ? !(element(at) == "H") : true
    # Check if the gmx executable exists
    if isnothing(Sys.which(gmx))
        throw(ArgumentError("""\n
            Could not find GROMACS `$gmx` executable. Add the `gmx` command to the path or provide it explicitly with the `gmx` keyword argument.

        """))
    end
    index_file = tempname() * ".ndx"
    sasa_file = tempname() * ".xvg"
    p = read_pdb(pdbname, at -> isprotein(at) && keepH(at))
    inds_protein = index.(p)
    inds_sidechain = index.(select(p, at -> resname(at) == rname && sidechain(at)))
    inds_backbone = index.(select(p, at -> resname(at) == rname && backbone(at)))
    open(index_file, "w") do io
        write(io, "[ SC ]\n")
        for i in inds_sidechain
            write(io, "$i\n")
        end
        write(io, "[ BB ]\n")
        for i in inds_backbone
            write(io, "$i\n")
        end
        write(io, "[ PROT ]\n")
        for i in inds_protein
            write(io, "$i\n")
        end
    end
    if length(inds_sidechain) == 0 # special case for GLY
        sasa_sc = 0.0
        try
            run(pipeline(`$gmx sasa -s $pdbname -probe 0.14 -ndots $n_dots -surface PROT -output BB -n $index_file -o $sasa_file`, stdout=devnull, stderr=devnull))
        catch
            error("Error running $gmx sasa for of $resname in $pdbname")
        end
        sasa_bb = read_gmx_delta_sasa_per_restype_values(sasa_file, 1)[1]
    else
        try
            run(pipeline(`$gmx sasa -s $pdbname -probe 0.14 -ndots $n_dots -surface PROT -output SC BB -n $index_file -o $sasa_file`, stdout=devnull, stderr=devnull))
        catch
            error("Error running $gmx sasa for of $resname in $pdbname")
        end
        sasa_temp = read_gmx_delta_sasa_per_restype_values(sasa_file, 2)
        sasa_sc, sasa_bb = sasa_temp[1], sasa_temp[2]
    end
    rm(sasa_file, force=true)
    rm(index_file, force=true)
    return sasa_bb, sasa_sc
end

@testitem "gmx_delta_sasa_per_restype (deprecated)" begin
    # These tests run only if GROMACS (`gmx`) is available: it is not installed in CI.
    using LazyArtifacts
    using PDBTools: mvalue_delta_sasa, gmx_delta_sasa_per_restype
    gmx = Sys.which("gmx")
    if isnothing(gmx)
        @warn "gmx executable not available: gmx_delta_sasa_per_restype tests won't be run"
    else
        @test_throws ArgumentError gmx_delta_sasa_per_restype(;
            native_pdb=PDBTools.test_file("1MJC_native.pdb"),
            desnat_pdb=PDBTools.test_file("1MJC_straight.pdb"),
            gmx="/tmp/gmx_fake",
        )

        MJC_native = read_pdb(PDBTools.test_file("1MJC_native.pdb"), "protein")
        sasa_1MJC = gmx_delta_sasa_per_restype(;
            native_pdb=PDBTools.test_file("1MJC_native.pdb"),
            desnat_pdb=PDBTools.test_file("1MJC_straight.pdb"),
        )
        r_1MJC = mvalue_delta_sasa(; atoms=MJC_native, sasas=sasa_1MJC)
        @test isapprox(r_1MJC.tot, -1.2335; rtol=1e-1)
        @test isapprox(r_1MJC.bb, -0.7333; rtol=1e-1)
        @test isapprox(r_1MJC.sc, -0.50025; rtol=1e-1)
        sasa_1MJC = gmx_delta_sasa_per_restype(;
            native_pdb=PDBTools.test_file("1MJC_native.pdb"),
            desnat_pdb=PDBTools.test_file("1MJC_straight.pdb"),
            ignore_hydrogen=false,
        )
        r_1MJC = mvalue_delta_sasa(; atoms=MJC_native, sasas=sasa_1MJC)
        @test isapprox(r_1MJC.tot, -0.98749; rtol=1e-1)
        @test isapprox(r_1MJC.bb, -0.3014; rtol=1e-1)
        @test isapprox(r_1MJC.sc, -0.67707; rtol=1e-1)

        RN2_native = read_pdb(PDBTools.test_file("2RN2_native.pdb"), "protein")
        sasa_2RN2 = gmx_delta_sasa_per_restype(;
            native_pdb=PDBTools.test_file("2RN2_native.pdb"),
            desnat_pdb=PDBTools.test_file("2RN2_straight.pdb"),
            gmx=gmx,
        )
        @test length(sasa_2RN2) == 20
        r_2RN2 = mvalue_delta_sasa(; atoms=RN2_native, sasas=sasa_2RN2)
        @test isapprox(r_2RN2.tot, -3.17; rtol=1e-1)
        @test isapprox(r_2RN2.bb, -1.94; rtol=1e-1)
        @test isapprox(r_2RN2.sc, -1.23; rtol=1e-1)
        sasa_2RN2 = gmx_delta_sasa_per_restype(;
            native_pdb=PDBTools.test_file("2RN2_native.pdb"),
            desnat_pdb=PDBTools.test_file("2RN2_straight.pdb"),
            ignore_hydrogen=false,
        )
        r_2RN2 = mvalue_delta_sasa(; atoms=RN2_native, sasas=sasa_2RN2)
        @test isapprox(r_2RN2.tot, -2.54; rtol=1e-1)
        @test isapprox(r_2RN2.bb, -0.88; rtol=1e-1)
        @test isapprox(r_2RN2.sc, -1.66; rtol=1e-1)
    end
end
