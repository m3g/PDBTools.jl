const creamer_atomic_radii = Dict{StringType, Float32}(
    "Nsp2" => 1.64,
    "Nsp3" => 1.64,
    "Csp2" => 1.76,
    "Csp3" => 1.88,
    "Osp2" => 1.42,
    "Osp3" => 1.46,
    "Ssp3" => 1.46,
    "H" => 0.0,
)

const PDB_ATOM_HYBRIDIZATION = Dict{StringType, Dict{StringType,StringType}}(
    # --- Backbone Atoms (Common to all amino acids, except Gly) ---
    "bb" => Dict(
        "N"    => "Nsp2", # Backbone Amide Nitrogen (peptide bond is N-C=O, typically sp2/sp3 intermediate)
        "CA"   => "Csp3", # Alpha Carbon
        "C"    => "Csp2", # Backbone Carbonyl Carbon
        "O"    => "Osp2", # Backbone Carbonyl Oxygen
        "OXT"  => "Osp2", # C-terminal Oxygen (in the last residue)
        "OT1"  => "Osp2", # C-terminal Oxygen (in the last residue)
        "OT2"  => "Osp2", # C-terminal Oxygen (in the last residue)
    ),

    # --- Glycine (GLY) Side Chain ---
    # Glycine has no side chain beyond CA, but sometimes the H is named H2 or H3
    "GLY" => Dict(
        "HA2"  => "H",    # Alpha Hydrogen
        "HA3"  => "H",    # Alpha Hydrogen
    ),

    # --- Alanine (ALA) Side Chain ---
    "ALA" => Dict(
        "CB"   => "Csp3", # Beta Carbon
    ),
 
    # --- Valine (VAL) Side Chain ---
    "VAL" => Dict(
        "CB"   => "Csp3",
        "CG1"  => "Csp3",
        "CG2"  => "Csp3",
    ),

    # --- Leucine (LEU) Side Chain ---
    "LEU" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp3",
        "CD1"  => "Csp3",
        "CD2"  => "Csp3",
    ),

    # --- Isoleucine (ILE) Side Chain ---
    "ILE" => Dict(
        "CB"   => "Csp3",
        "CG1"  => "Csp3",
        "CG2"  => "Csp3",
        "CD"   => "Csp3",
        "CD1"  => "Csp3",
    ),

    # --- Serine (SER) Side Chain ---
    "SER" => Dict(
        "CB"   => "Csp3",
        "OG"   => "Osp3", # Gamma Oxygen (Alcohol)
    ),

    # --- Threonine (THR) Side Chain ---
    "THR" => Dict(
        "CB"   => "Csp3",
        "OG1"  => "Osp3", # Gamma 1 Oxygen (Alcohol)
        "CG2"  => "Csp3",
    ),

    # --- Cysteine (CYS) Side Chain ---
    "CYS" => Dict(
        "CB"   => "Csp3",
        "SG"   => "Ssp3", # Gamma Sulfur (Thiol)
    ),

    # --- Methionine (MET) Side Chain ---
    "MET" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp3",
        "SD"   => "Ssp3", # Delta Sulfur (Thioether)
        "CE"   => "Csp3",
    ),

    # --- Proline (PRO) Side Chain (Cyclic with backbone N) ---
    # The backbone N in Proline is an sp3 tertiary amine.
    # N is often considered Nsp3 in Proline, unlike other residues.
    "PRO" => Dict(
        "N"   => "Nsp3", # Specific for Proline's N
        "CB"  => "Csp3",
        "CG"  => "Csp3",
        "CD"  => "Csp3", # Delta Carbon
    ),

    # --- Phenylalanine (PHE) Side Chain ---
    "PHE" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp2", # Benzene ring attachment
        "CD1"  => "Csp2", # Ring Carbons
        "CD2"  => "Csp2",
        "CE1"  => "Csp2",
        "CE2"  => "Csp2",
        "CZ"   => "Csp2",
    ),

    # --- Tyrosine (TYR) Side Chain ---
    "TYR" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp2",
        "CD1"  => "Csp2",
        "CD2"  => "Csp2",
        "CE1"  => "Csp2",
        "CE2"  => "Csp2",
        "CZ"   => "Csp2",
        "OH"   => "Osp3", # Hydroxyl Oxygen
    ),

    # --- Tryptophan (TRP) Side Chain (Indole Ring) ---
    "TRP" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp2",
        "CD1"  => "Csp2",
        "CD2"  => "Csp2",
        "NE1"  => "Nsp2", # Indole Nitrogen
        "CE2"  => "Csp2",
        "CE3"  => "Csp2",
        "CZ2"  => "Csp2",
        "CZ3"  => "Csp2",
        "CH2"  => "Csp2",
    ),

    # --- Aspartic Acid (ASP) Side Chain ---
    "ASP" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp2", # Carboxyl Carbon
        "OD1"  => "Osp2", # Carboxyl Oxygen
        "OD2"  => "Osp2", # Carboxyl Oxygen
    ),

    # --- Asparagine (ASN) Side Chain ---
    "ASN" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp2", # Amide Carbon
        "OD1"  => "Osp2", # Amide Oxygen
        "ND2"  => "Nsp3", # Amide Nitrogen (can be Nsp2 in some contexts, but Nsp3 is better for simple type)
    ),

    # --- Glutamic Acid (GLU) Side Chain ---
    "GLU" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp3",
        "CD"   => "Csp2", # Carboxyl Carbon
        "OE1"  => "Osp2", # Carboxyl Oxygen
        "OE2"  => "Osp2", # Carboxyl Oxygen
    ),

    # --- Glutamine (GLN) Side Chain ---
    "GLN" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp3",
        "CD"   => "Csp2", # Amide Carbon
        "OE1"  => "Osp2", # Amide Oxygen
        "NE2"  => "Nsp3", # Amide Nitrogen (can be Nsp2 in some contexts, but Nsp3 is better for simple type)
    ),

    # --- Lysine (LYS) Side Chain (Like your example) ---
    "LYS" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp3",
        "CD"   => "Csp3",
        "CE"   => "Csp3",
        "NZ"   => "Nsp3", # Zeta Nitrogen (Primary Amine/Ammonium)
    ),

    # --- Arginine (ARG) Side Chain (Guanidinium Group) ---
    "ARG" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp3",
        "CD"   => "Csp3",
        "NE"   => "Nsp2", # Epsilon Nitrogen (Part of the delocalized system)
        "CZ"   => "Csp2", # Guanidinium Carbon
        "NH1"  => "Nsp2", # Guanidinium Nitrogen
        "NH2"  => "Nsp2", # Guanidinium Nitrogen
    ),

    # --- Histidine (HIS) Side Chain (Imidazole Ring) ---
    "HIS" => Dict(
        "CB"   => "Csp3",
        "CG"   => "Csp2", # Ring attachment
        "ND1"  => "Nsp2", # Delta 1 Nitrogen (Pyridine-like, or N-H Imidazole form)
        "CE1"  => "Csp2",
        "NE2"  => "Nsp2", # Epsilon 2 Nitrogen (Pyrole-like, or N Imidazole form)
        "CD2"  => "Csp2",
    ),
)

function creamer_atom_type(at::Atom)
    rname = threeletter(resname(at))
    atname = name(at)
    if haskey(PDB_ATOM_HYBRIDIZATION, rname)
        if haskey(PDB_ATOM_HYBRIDIZATION[rname], atname)
            return PDB_ATOM_HYBRIDIZATION[rname][atname]
        elseif haskey(PDB_ATOM_HYBRIDIZATION["bb"], atname)
            return PDB_ATOM_HYBRIDIZATION["bb"][atname]
        elseif element(at) == "H" 
            return "H"
        end
    end
    throw(ArgumentError(("""\n
        Could not determine Creamer united atom type for $atname of residue $(resname(at))

    """)))
end

const creamer_sasas = Dict{
    String,
    @NamedTuple{bb_lower::Float32, bb_upper::Float32, sc_lower::Float32, sc_upper::Float32}
}(
    "ALA" => (bb_lower=19.8, bb_upper=35.9, sc_lower=46.6, sc_upper=63.6),
    "ARG" => (bb_lower=17.1, bb_upper=33.0, sc_lower=156.9, sc_upper=185.3),
    "ASN" => (bb_lower=17.6, bb_upper=32.7, sc_lower=84.5, sc_upper=95.6),
    "ASP" => (bb_lower=18.1, bb_upper=33.9, sc_lower=79.2, sc_upper=94.8),
    "CYS" => (bb_lower=18.2, bb_upper=34.5, sc_lower=62.9, sc_upper=83.0),
    "GLN" => (bb_lower=17.2, bb_upper=33.4, sc_lower=105.0, sc_upper=128.7),
    "GLU" => (bb_lower=17.9, bb_upper=33.5, sc_lower=102.8, sc_upper=123.9),
    "GLY" => (bb_lower=54.6, bb_upper=75.7, sc_lower=0.0, sc_upper=0.0),
    "HIS" => (bb_lower=14.9, bb_upper=33.4, sc_lower=103.9, sc_upper=119.1),
    "ILE" => (bb_lower=15.2, bb_upper=24.7, sc_lower=100.1, sc_upper=134.1),
    "LEU" => (bb_lower=14.7, bb_upper=30.7, sc_lower=101.4, sc_upper=117.7),
    "LYS" => (bb_lower=18.3, bb_upper=33.8, sc_lower=142.5, sc_upper=158.8),
    "MET" => (bb_lower=16.7, bb_upper=33.8, sc_lower=105.3, sc_upper=139.5),
    "PHE" => (bb_lower=15.3, bb_upper=33.3, sc_lower=118.7, sc_upper=139.8),
    "PRO" => (bb_lower=18.9, bb_upper=26.1, sc_lower=83.5, sc_upper=90.5),
    "SER" => (bb_lower=23.8, bb_upper=35.0, sc_lower=59.7, sc_upper=73.3),
    "THR" => (bb_lower=18.6, bb_upper=29.5, sc_lower=77.3, sc_upper=91.2),
    "TRP" => (bb_lower=15.1, bb_upper=32.0, sc_lower=154.7, sc_upper=158.4),
    "TYR" => (bb_lower=17.7, bb_upper=33.5, sc_lower=131.0, sc_upper=152.3),
    "VAL" => (bb_lower=15.9, bb_upper=24.9, sc_lower=81.8, sc_upper=110.9),
)

struct _Selector{F,T} <: Function
    f::F
    residue::Ref{T}
end
(s::_Selector)(at::Atom) = s.f(at) && (at in s.residue[])

"""
    creamer_delta_sasa(atoms::AbstractVector{<:Atom})

Computes, for a vector of protein atoms, the predicted changes in SASA upon denaturation, using the 
Creamer model. Returns a dictionary that can be directly used as input to the `mvalue` function.

Reference:

Creamer TP, Srinivasan R, Rose GD. Modeling unfolded states of proteins and peptides. II. Backbone solvent accessibility. 
*Biochemistry.* 1997;36:2832–2835. doi: 10.1021/bi962819o.

# Example:

```julia
using PDBTools
using PDBTools: mvalue_delta_sasa, creamer_delta_sasa
ats = read_pdb("protein.pdb")
m = mvalue_delta_sasa(;
    model=AutonBolen,
    cosolvent="urea",
    atoms=ats,
    sasas=creamer_delta_sasa(ats)
)
```

The output is a tuple where `m.tot`, `m.bb`, `m.sc` are the total, backbone, and sidechain contributions for the denaturation
transfer free energy, and `m.restype` is a dictionary where each key contains the backbone and sidechain contributions 
for each residue type.

"""
function creamer_delta_sasa(atoms::AbstractVector{<:Atom})
    sasas = Dict{String,Dict}()
    sasa_atoms = sasa_particles(atoms;
        atom_type = creamer_atom_type,
        atom_radius_from_type = type -> creamer_atomic_radii[type]
    )
    sel_bb = _Selector(isbackbone, Ref(first(eachresidue(atoms))))
    sel_sc = _Selector(issidechain, Ref(first(eachresidue(atoms))))
    for res in eachresidue(atoms)
        rname = threeletter(resname(res))
        if !haskey(sasas, rname)
            sasas[rname] = Dict(:sc => (0.0, 0.0, 0.0), :bb => (0.0, 0.0, 0.0))
        end
    end
    for res in eachresidue(atoms)
        sel_bb.residue[] = res
        sel_sc.residue[] = res
        sasa_res_bb = sasa(sasa_atoms, sel_bb)
        sasa_res_sc = sasa(sasa_atoms, sel_sc)
        rname = threeletter(resname(res))
        cr = creamer_sasas[rname]
        csc = sasas[rname][:sc]
        cbb = sasas[rname][:bb]
        sasas[rname] = Dict(
            :sc => (
                csc[1] + cr.sc_lower - sasa_res_sc,
                csc[2] + 0.5 * (cr.sc_lower + cr.sc_upper) - sasa_res_sc,
                csc[3] + cr.sc_upper - sasa_res_sc,
            ),
            :bb => (
                cbb[1] + cr.bb_lower - sasa_res_bb,
                cbb[2] + 0.5 * (cr.bb_lower + cr.bb_upper) - sasa_res_bb,
                cbb[3] + cr.bb_upper - sasa_res_bb,
            )
        )
    end
    return sasas
end
