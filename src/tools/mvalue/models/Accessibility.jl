#=

Accessibility model, where:

tfe_sc = (1/f_sc) * (GTFEapp - γGLY + (1 + f_bb)*tfe_bb)

=#

export Accessibility
struct Accessibility <: MValueModel end

# Do not user underscores (_) in the following names:
const cosolvent_column_Accessibility = OrderedDict(
    "tmao" => 1,
    "sarcosine" => 2,
    "betaine" => 3,
    "proline" => 4,
    "sorbitol" => 5,
    "sucrose" => 6,
    "urea" => 7,
    "glycerol" => 8,
    "trehalose" => 9,
)

#
# ASA of backbones and side-chains computed: 
# bb: accessible surface area of a backbone of the center peptide in a 
#     tri-peptide, disconsidering the side-chains of the first and last residues 
#     emulating GXG.
# bb_pure: Accessible surface area of the central backbone of XYZ, ignoring
#          all backbones.
# f_bb = bb / bb_pure - the fraction of the BB protected by the central side chain.
#
# sc: accessible surface area of the side-chain of an isolated residue.
# sc_pure: accessible surface are of the side-chain detached from the backbone.
# f_sc = sc / sc_pure: the fraction of the SC that is protected shielded by the BB.
#
# n: number of entries in the dataase: (XYZ) fragments with a given type of Y
#
# Database: All CATH 20 domains, as downloaded in 26/05/2026
#
const f_acc = OrderedDict{String, OrderedDict{String, Float32}}(
  "ALA" => OrderedDict("n"=>178636, "sc"=>62.147, "sc_pure"=>120.728, "bb"=>52.93, "bb_pure"=>87.459, "f_bb"=>0.605198, "f_sc"=>0.514768),
  "PHE" => OrderedDict("n"=>91160, "sc"=>180.799, "sc_pure"=>237.736, "bb"=>42.2336, "bb_pure"=>86.2234, "f_bb"=>0.489816, "f_sc"=>0.760502),
  "LEU" => OrderedDict("n"=>212344, "sc"=>144.686, "sc_pure"=>200.589, "bb"=>42.3231, "bb_pure"=>86.9463, "f_bb"=>0.486774, "f_sc"=>0.721302),
  "ILE" => OrderedDict("n"=>126302, "sc"=>144.897, "sc_pure"=>201.34, "bb"=>40.0158, "bb_pure"=>86.3755, "f_bb"=>0.463277, "f_sc"=>0.719663),
  "VAL" => OrderedDict("n"=>152561, "sc"=>120.343, "sc_pure"=>177.027, "bb"=>41.3333, "bb_pure"=>86.1337, "f_bb"=>0.479873, "f_sc"=>0.6798),
  "PRO" => OrderedDict("n"=>97658, "sc"=>115.728, "sc_pure"=>175.361, "bb"=>42.857, "bb_pure"=>90.3427, "f_bb"=>0.474383, "f_sc"=>0.659941),
  "MET" => OrderedDict("n"=>36817, "sc"=>155.679, "sc_pure"=>212.483, "bb"=>44.0675, "bb_pure"=>86.7157, "f_bb"=>0.508184, "f_sc"=>0.732663),
  "TRP" => OrderedDict("n"=>31439, "sc"=>222.825, "sc_pure"=>279.907, "bb"=>41.012, "bb_pure"=>86.8561, "f_bb"=>0.472184, "f_sc"=>0.796068),
  "GLY" => OrderedDict("n"=>151206, "sc"=>0.0, "sc_pure"=>0.0, "bb"=>86.3549, "bb_pure"=>86.3549, "f_bb"=>1.0, "f_sc"=>1.0),
  "SER" => OrderedDict("n"=>129548, "sc"=>83.0034, "sc_pure"=>141.191, "bb"=>48.602, "bb_pure"=>86.5017, "f_bb"=>0.561862, "f_sc"=>0.587879),
  "THR" => OrderedDict("n"=>118271, "sc"=>111.981, "sc_pure"=>168.956, "bb"=>42.4243, "bb_pure"=>85.8036, "f_bb"=>0.494435, "f_sc"=>0.662778),
  "TYR" => OrderedDict("n"=>78108, "sc"=>198.268, "sc_pure"=>255.194, "bb"=>42.4376, "bb_pure"=>86.1911, "f_bb"=>0.492366, "f_sc"=>0.77693),
  "GLN" => OrderedDict("n"=>82374, "sc"=>151.399, "sc_pure"=>208.416, "bb"=>44.3666, "bb_pure"=>87.068, "f_bb"=>0.509563, "f_sc"=>0.726426),
  "ASN" => OrderedDict("n"=>93050, "sc"=>126.426, "sc_pure"=>184.088, "bb"=>44.2684, "bb_pure"=>87.562, "f_bb"=>0.505566, "f_sc"=>0.686771),
  "ASP" => OrderedDict("n"=>125402, "sc"=>124.667, "sc_pure"=>182.223, "bb"=>44.11, "bb_pure"=>87.8884, "f_bb"=>0.501887, "f_sc"=>0.684146),
  "GLU" => OrderedDict("n"=>141726, "sc"=>149.409, "sc_pure"=>206.402, "bb"=>44.8279, "bb_pure"=>87.6843, "f_bb"=>0.511242, "f_sc"=>0.723873),
  "HIS" => OrderedDict("n"=>50912, "sc"=>157.599, "sc_pure"=>214.954, "bb"=>43.8755, "bb_pure"=>86.5002, "f_bb"=>0.50723, "f_sc"=>0.733173),
  "LYS" => OrderedDict("n"=>117097, "sc"=>165.952, "sc_pure"=>223.125, "bb"=>46.1709, "bb_pure"=>87.3999, "f_bb"=>0.528271, "f_sc"=>0.743764),
  "ARG" => OrderedDict("n"=>109623, "sc"=>198.673, "sc_pure"=>256.034, "bb"=>45.7085, "bb_pure"=>87.0432, "f_bb"=>0.525124, "f_sc"=>0.775961),
  "CYS" => OrderedDict("n"=>29486, "sc"=>103.039, "sc_pure"=>160.6, "bb"=>46.205, "bb_pure"=>85.928, "f_bb"=>0.537717, "f_sc"=>0.641584),
)

#
# Backbone accessibility parameter - [-Infty,1]
#
# if acc == 1 -> the accessibility is 1.0 (full access)
# if acc == 0 -> the accessibility is defined by the ratio of ASAs of actual and pure BBs
# if acc < 0 -> the accessibility is greater than that of the ASAs ratio
#
# acc(f_bb, x) = f_bb^(1 - x)
#
const acc = Dict{String,Float32}(
    "tmao" => 0.0,
    "sarcosine" => 0.0,
    "betaine" => 0.0,
    "proline" => 0.0,
    "glycerol" => 0.0,
    "sorbitol" => 0.0,
    "sucrose" => 0.0,
    "trehalose" => 0.0,
    "urea" => 1.0,
)

function model_combination_rule(::Type{Accessibility}, cosolvent, restype)
    col = data_col[cosolvent_column_Accessibility[cosolvent]]
    # united model: all bb ASA contributions are the same
    bb_contribution = GTFEapp["BB"][col] / f_acc["GLY"]["bb_pure"]
    sc_contribution = if restype == "GLY"
        0.0f0
    else
        f_sc = f_acc[restype]["f_sc"] # accessibility of the side-chain
        f_bb = f_acc[restype]["f_bb"] # accessibility of the backbone
        _acc = f_bb^(1 - acc[cosolvent]) # f_bb (protectant) or 1 (urea)
        GTFE_sc = GTFEapp[restype][col] # apparent free energy of the backbone
        tfe_bb = GTFEapp["BB"][col] # apparent free energy of the backbone
        γAA = γ[restype][col] # amino-acid activity correction (usually ignored)
        γG = γ["GLY"][col] # gly-activitity correction - affects all side-chains
        TFE_sc = (1/f_sc) * (GTFE_sc + γAA - γG + tfe_bb * (1 - _acc))
        TFE_sc / f_acc[restype]["sc_pure"]
    end
    return bb_contribution, sc_contribution
end