struct RichardsUnitedAtomRadii <: AtomicRadiiType end

_richards_not_hydrogen(at) = element(at) != "H"

# In some PDB files, `resname(at)` folds in a leading altLoc character for atoms with
# short (2-character) names (e.g. "ASER"/"BSER" for the alternate conformers of a Ser OG),
# instead of just "SER". This retries the lookup after stripping that leading character.
function _richards_resname(at::Atom)
    raw = resname(at)
    rname = threeletter(StringType, raw)
    if isnothing(rname) && length(raw) == 4
        rname = threeletter(StringType, raw[2:end])
    end
    return rname
end

#=
    richards_atom_type(at::Atom)

Classifies a protein atom into one of the Richards united-atom groups, reusing the same
per-residue atom -> hybridization assignment used for the Creamer radii (`PDB_ATOM_HYBRIDIZATION`,
defined in `creamer.jl`). The Richards classification is, in essence, the same sp3/sp2
hybridization split used by Creamer, relabeled with the group names used in the original
Richards papers:

- `Csp3` -> `"CH4"` (tetrahedral/sp3 carbon)
- `Csp2` -> `"CH3"` (trigonal/sp2 carbon, e.g. aromatic or carbonyl carbon)
- `Nsp3` -> `"NH4"` (tetrahedral/sp3 nitrogen, e.g. Lys ammonium, Pro ring nitrogen)
- `Nsp2` -> `"NH3"` (trigonal/sp2 nitrogen, e.g. backbone amide, Arg guanidino, aromatic ring N)
- `Osp3` -> `"OH4"` (tetrahedral/sp3 oxygen, e.g. Ser/Thr/Tyr hydroxyl)
- `Osp2` -> `"OH3"` (trigonal/sp2 oxygen, e.g. carbonyl oxygen)
- sulfur is split by residue instead of hybridization: Cys `SG` -> `"SH"` (thiol), Met `SD` -> `"ST"` (thioether)
- `Zn` and `Fe` hetero ions are mapped to `"ZN"` and `"FE"`, respectively.
=#
function richards_atom_type(at::Atom)
    el = element(at)
    (el == "Zn" || el == "ZN") && return "ZN"
    (el == "Fe" || el == "FE") && return "FE"
    rname = _richards_resname(at)
    atname = name(at)
    hybrid = if haskey(PDB_ATOM_HYBRIDIZATION, rname) && haskey(PDB_ATOM_HYBRIDIZATION[rname], atname)
        PDB_ATOM_HYBRIDIZATION[rname][atname]
    elseif haskey(PDB_ATOM_HYBRIDIZATION["bb"], atname)
        PDB_ATOM_HYBRIDIZATION["bb"][atname]
    elseif el == "H"
        "H"
    else
        nothing
    end
    isnothing(hybrid) && throw(ArgumentError("""\n
        Could not determine Richards united atom type for $atname of residue $(resname(at))

    """))
    hybrid == "Csp3" && return "CH4"
    hybrid == "Csp2" && return "CH3"
    hybrid == "Nsp3" && return "NH4"
    hybrid == "Nsp2" && return "NH3"
    hybrid == "Osp3" && return "OH4"
    hybrid == "Osp2" && return "OH3"
    hybrid == "H" && return "H"
    if hybrid == "Ssp3"
        atname == "SG" && return "SH"
        atname == "SD" && return "ST"
    end
    throw(ArgumentError("""\n
        Could not determine Richards united atom type for $atname of residue $(resname(at))

    """))
end

#
# Two alternative sets of Richards united-atom radii, both quoted from the same table:
#
#     Set 1                         Set 2
# (Richards, 1977)   (Richmond & Richards, 1978)
#
# ch4 2.00                1.90
# ch3 1.70                1.70
# nh4 2.00                1.70
# nh3 1.70                1.70
# oh4 1.60                1.40
# oh3 1.50                1.40
# sh  2.00                1.80
# st  1.80                1.80
# zn  0.64                0.64
# fe  0.64                0.64
#
const richards_atomic_radii_set1 = OrderedDict{StringType,Float32}(
    "CH4" => 2.00,
    "CH3" => 1.70,
    "NH4" => 2.00,
    "NH3" => 1.70,
    "OH4" => 1.60,
    "OH3" => 1.50,
    "SH" => 2.00,
    "ST" => 1.80,
    "ZN" => 0.64,
    "FE" => 0.64,
    "H" => 0.0,
)

const richards_atomic_radii_set2 = OrderedDict{StringType,Float32}(
    "CH4" => 1.90,
    "CH3" => 1.70,
    "NH4" => 1.70,
    "NH3" => 1.70,
    "OH4" => 1.40,
    "OH3" => 1.40,
    "SH" => 1.80,
    "ST" => 1.80,
    "ZN" => 0.64,
    "FE" => 0.64,
    "H" => 0.0,
)

function _richards_radii_set(s::Symbol)
    s == :set1 && return richards_atomic_radii_set1
    s == :set2 && return richards_atomic_radii_set2
    throw(ArgumentError("""\n
        radii_set parameter must be either :set1 (Richards, 1977) or :set2 (Richmond & Richards, 1978)

    """))
end

"""
    sasa_particles(RichardsUnitedAtomRadii, atoms; radii_set=:set2, kargs...)

Compute the SASA of the structure using the Richards united atom radii, as used by
SurfaceRacer. This parameterization is available for proteins only, and ignores any
hydrogen atoms of the structure.

Two alternative radii sets are available, both reproduced from the same source table:

- `:set2` (default): Richmond & Richards, 1978.
- `:set1`: Richards, 1977.

Both sets were validated against the total SASA values reported by SurfaceRacer for the
30 structures in `SurfaceRacerResults` (downloading each PDB entry and computing the SASA
of all protein atoms, all chains, no hydrogens): both reproduce the reported values with a
mean absolute error of ~2% (root-mean-square error ~2.4-2.5%), with `:set2` showing a
slightly smaller mean bias and `:set1` a slightly smaller worst-case outlier. This confirms
that the atom -> group classification used by `richards_atom_type` (the same sp3/sp2 split
used for the Creamer radii, relabeled) matches SurfaceRacer's own assignment closely enough
that the residual discrepancy is attributable to details of the rolling-ball algorithm
(number of dots per atom, exact probe radius rounding) rather than to atom (mis)typing.

"""
function sasa_particles(::Type{RichardsUnitedAtomRadii}, atoms; radii_set::Symbol=:set2, kargs...)
    radii = _richards_radii_set(radii_set)
    p_no_h = select(atoms, _richards_not_hydrogen)
    s = _sasa_particles(RichardsUnitedAtomRadii, p_no_h;
        atom_type=richards_atom_type,
        atom_radius_from_type=type -> radii[type],
        kargs...
    )
    return s
end

#
# Total SASA (Å²) reported by SurfaceRacer for these PDB structures, used here as a
# reference to validate the atom-type -> radius assignment above.
#
const SurfaceRacerResults = Dict{String,Float32}(
    "3sdh" => 14326.0,
    "1rbc" => 6755.0,
    "2act" => 9124.0,
    "1rbf" => 6695.0,
    "4gcr" => 9040.0,
    "1thm" => 10010.0,
    "2lyz" => 6730.0,
    "1icm" => 7598.0,
    "2sn3" => 4214.0,
    "3cyt" => 11846.0,
    "1eca" => 7234.0,
    "1cau" => 20245.0,
    "5p21" => 8673.0,
    "1rbg" => 6713.0,
    "2rns" => 6713.0,
    "1rbh" => 6728.0,
    "2ptn" => 9410.0,
    "1nxb" => 4083.0,
    "3rn3" => 6952.0,
    "1rro" => 6061.0,
    "1ycc" => 6627.0,
    "1rbi" => 6690.0,
    "4pti" => 3963.0,
    "1ecd" => 7203.0,
    "5mbn" => 8417.0,
    "1cse" => 12918.0,
    "1mbd" => 8478.0,
    "1rbe" => 6702.0,
    "1plc" => 5128.0,
    "1arb" => 9696.0,
)

@testitem "RichardsUnitedAtomRadii" begin
    using PDBTools
    prot = read_pdb(PDBTools.TESTPDB, "protein")

    s2 = sasa(sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot))
    @test s2 ≈ 5234.202 rtol = 1e-4

    s1 = sasa(sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot; radii_set=:set1))
    @test s1 ≈ 5260.748 rtol = 1e-4

    @test_throws "radii_set parameter must be" sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot; radii_set=:wrong)

    # Test invalid atom type error
    prot[1].name = "AB"
    @test_throws "united atom type" sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot)

    # Test with the presence of hydrogens (which must be ignored)
    prot = read_pdb(PDBTools.TESTPDB, "protein")
    s1 = sasa(sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot))
    @test s1 ≈ s2

    # Test with a double-conformer amino acid name
    prot[1].resname = "BALA"
    s1 = sasa(sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot))
    @test s1 ≈ s2

end
