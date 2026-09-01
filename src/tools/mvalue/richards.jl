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
    if isnothing(hybrid) 
        throw(ArgumentError("""\n
            Could not determine Richards united atom type for $atname of residue $(resname(at))
    
        """))
    end
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
    # Unreachable with the current PDB_ATOM_HYBRIDIZATION table (every "Ssp3" entry is
    # ("CYS","SG") or ("MET","SD")); kept as a guard against a silent `nothing` return
    # if that table is ever extended with another Ssp3-mapped atom name.
    throw(ArgumentError("""\n
        Could not determine Richards united atom type for $atname of residue $(resname(at))

    """))
end

#
# Three alternative sets of Richards united-atom radii, corresponding to the three
# choices offered by SurfaceRacer itself ("1 - Richards (1977)" and "2 - Chothia
# (1976)" in its own prompt, plus the Richmond & Richards variant used elsewhere):
#
#     Set 1                Set 2                    Set 3
# (Richards, 1977)   (Richmond & Richards, 1978)   (Chothia, 1976)
#
# ch4 2.00                1.90                       1.87
# ch3 1.70                1.70                       1.76
# nh4 2.00                1.70                       1.50
# nh3 1.70                1.70                       1.65
# oh4 1.60                1.40                       1.40
# oh3 1.50                1.40                       1.40
# sh  2.00                1.80                       1.85
# st  1.80                1.80                       1.85
# zn  0.64                0.64                       0.64
# fe  0.64                0.64                       0.64
#
# Set 3 was validated directly against a local run of SurfaceRacer 5.0 (Tsodikov,
# Record & Sergeev, 2002) on the biological tetramer assembly of 3CNA (Hardman &
# Ainsworth, 1972), with its own "2 - Chothia (1976)" option and a 1.4 Å probe: total
# SASA 35477.35 Å² there vs. ~35480-35575 Å² here depending on `n_dots` (<0.3%
# difference, consistent with an exact analytical algorithm vs. this module's
# numerical, dot-based one). See the "RichardsUnitedAtomRadii validated against
# SurfaceRacer" testitem below.
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

# From Chothia 1976 (J. Mol. Biol. 105, 1-14)
const richards_atomic_radii_set3 = OrderedDict{StringType,Float32}(
    "CH4" => 1.87,
    "CH3" => 1.76,
    "NH4" => 1.50,
    "NH3" => 1.65,
    "OH4" => 1.40,
    "OH3" => 1.40,
    "SH" => 1.85,
    "ST" => 1.85,
    "ZN" => 0.64,
    "FE" => 0.64,
    "H" => 0.0,
)

function _richards_radii_set(s::Symbol)
    s == :set1 && return richards_atomic_radii_set1
    s == :set2 && return richards_atomic_radii_set2
    s == :set3 && return richards_atomic_radii_set3
    throw(ArgumentError("""\n
        radii_set parameter must be one of
            :set1 (Richards, 1977) 
            :set2 (Richmond & Richards, 1978)
            :set3 (Chothia, 1976)

    """))
end

"""
    sasa_particles(RichardsUnitedAtomRadii, atoms; radii_set=:set2, kargs...)

Compute the SASA of the structure using the Richards united atom radii, as used by
SurfaceRacer. This parameterization is available for proteins only, and ignores any
hydrogen atoms of the structure.

Three alternative radii sets are available, matching SurfaceRacer's own two built-in
choices plus a third variant used elsewhere in the transfer-model literature:

- `:set2` (default): Richmond & Richards, 1978.
- `:set1`: Richards, 1977 (SurfaceRacer's own "1 - Richards (1977)" option).
- `:set3`: Chothia, J. Mol. Biol., 1976 (SurfaceRacer's own "2 - Chothia (1976)" option).

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

    s3 = sasa(sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot; radii_set=:set3))
    @test s3 ≈ 5246.431 rtol = 1e-4

    @test_throws "radii_set parameter must be" sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot; radii_set=:wrong)

    # Test invalid atom type error
    prot[1].name = "AB"
    @test_throws "united atom type" sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot)
    prot[1].name = "SF"
    @test_throws "united atom type" sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot)

    # Test with the presence of hydrogens (which must be ignored)
    prot = read_pdb(PDBTools.TESTPDB, "protein")
    s1 = sasa(sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot))
    @test s1 ≈ s2
    @test PDBTools.richards_atom_type(prot[2]) == "H"

    # Test with a double-conformer amino acid name
    prot[1].resname = "BALA"
    s1 = sasa(sasa_particles(PDBTools.RichardsUnitedAtomRadii, prot))
    @test s1 ≈ s2

end

@testitem "RichardsUnitedAtomRadii validated against SurfaceRacer" begin
    using PDBTools

    # SurfaceRacer 5.0 (Tsodikov, Record & Sergeev, 2002), run locally with its own
    # "2 - Chothia (1976)" radii option and a 1.4 Å probe, on the biological tetramer
    # assembly of 3CNA (Hardman & Ainsworth, 1972; 4 chains, 7228 protein atoms):
    #
    #     tetramer:                        Total area = 35477.35 Å²
    #     each dimer (chain A+A-2 or A-3+A-4, related by an exact crystallographic
    #     2-fold, and therefore geometrically identical): Total area = 19954.90 Å² (both)
    #
    # `:set3` (with the corrected NH4 radius) reproduces the tetramer value to <0.3%,
    # validating the radii table and the dot-based SASA algorithm against the reference
    # program (residual difference expected: exact analytical vs. numerical dot-based
    # methods). Note `wget("3CNA")` only returns the asymmetric unit (a single chain);
    # the biological tetramer requires RCSB's separate assembly file, hence `assembly=1`.
    cna = wget("3CNA", "protein"; assembly=1)
    @test length(cna) == 7228
    @test Set(unique(chain.(cna))) == Set(["A", "A-2", "A-3", "A-4"])

    d1 = select(cna, at -> chain(at) in ("A", "A-2"))
    d2 = select(cna, at -> chain(at) in ("A-3", "A-4"))
    @test length(d1) == length(d2) == 3614

    s_tet = sasa(sasa_particles(RichardsUnitedAtomRadii, cna; radii_set=:set3))
    @test s_tet ≈ 35477.35 rtol = 0.005

    # Higher n_dots here only to demonstrate the exact crystallographic symmetry
    # (SASA(d1) == SASA(d2)) without dot-pattern-orientation noise; the package
    # default (n_dots=512) already agrees with SurfaceRacer to within ~0.2%.
    s_d1 = sasa(sasa_particles(RichardsUnitedAtomRadii, d1; radii_set=:set3, n_dots=5000))
    s_d2 = sasa(sasa_particles(RichardsUnitedAtomRadii, d2; radii_set=:set3, n_dots=5000))
    @test s_d1 ≈ s_d2 rtol = 0.002
    @test s_d1 ≈ 19954.90 rtol = 0.02
    @test s_d2 ≈ 19954.90 rtol = 0.02

    # Regression pin (default radii_set=:set2) for the tetramer -> 2-dimer dissociation
    # transfer free energies, computed exactly as Table 2 of Knowles et al. 2015
    # (10.1021/acs.biochem.5b00246) defines the "no conformational changes" ΔASA:
    # Δtfe = tfe(d1) + tfe(d2) - tfe(tetramer). These do NOT reproduce that paper's own
    # predicted m-values (traced to structure-preparation details not fully identified);
    # they are pinned here only to catch future regressions in this codebase.
    targets = Dict(
        "tetraeg" => 0.19622135,
        "urea" => -0.45273495,
        "glycerol" => 0.13829088,
        "proline" => 0.5135155,
        "betaine" => 0.469985,
    )
    for (s, target) in targets
        delta = transfer_free_energy(d1, s; model=MTRecord).tot +
                transfer_free_energy(d2, s; model=MTRecord).tot -
                transfer_free_energy(cna, s; model=MTRecord).tot
        @test delta ≈ target rtol = 1e-4
    end
end
