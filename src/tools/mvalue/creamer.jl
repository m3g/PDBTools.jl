
const creamer_atomic_radii = Dict{String, Float32}(
  "C_BB" => 1.79547,
  "O_BB" => 1.39941,
  "S_SC" => 1.73841,
  "O_SC" => 1.38058,
  "N_BB" => 1.76601,
  "N_SC" => 1.38848,
  "C_SC" => 1.88785,
)

function creamer_atom_type(at)
    name(at) == "N" && return "N_BB"
    name(at) == "CA" && return "C_BB"
    name(at) == "C" && return "C_BB"
    name(at) == "O" && return "O_BB"
    element(at) == "C" && return "C_SC"
    element(at) == "N" && return "N_SC"
    element(at) == "O" && return "O_SC"
    element(at) == "S" && return "S_SC"
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
        rname = resname(res)
        if !haskey(sasas, rname)
            sasas[rname] = Dict(:sc => (0.0, 0.0, 0.0), :bb => (0.0, 0.0, 0.0))
        end
    end
    for res in eachresidue(atoms)
        sel_bb.residue[] = res
        sel_sc.residue[] = res
        sasa_res_bb = sasa(sasa_atoms, sel_bb)
        sasa_res_sc = sasa(sasa_atoms, sel_sc)
        rname = resname(res)
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
