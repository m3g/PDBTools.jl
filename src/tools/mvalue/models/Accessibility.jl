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
    "ALA" => OrderedDict("n"=>179233, "sc"=>70.878, "sc_pure"=>135.348, "bb"=>49.0121, "bb_pure"=>88.8883, "f_bb"=>0.55139, "f_sc"=>0.523674),
    "PHE" => OrderedDict("n"=>91323, "sc"=>187.587, "sc_pure"=>248.851, "bb"=>38.6632, "bb_pure"=>87.374, "f_bb"=>0.442502, "f_sc"=>0.753811),
    "LEU" => OrderedDict("n"=>213010, "sc"=>159.14, "sc_pure"=>219.486, "bb"=>37.3229, "bb_pure"=>88.2639, "f_bb"=>0.422855, "f_sc"=>0.725056),
    "ILE" => OrderedDict("n"=>126589, "sc"=>159.817, "sc_pure"=>220.105, "bb"=>35.3226, "bb_pure"=>87.5728, "f_bb"=>0.403351, "f_sc"=>0.726092),
    "VAL" => OrderedDict("n"=>153028, "sc"=>133.941, "sc_pure"=>194.604, "bb"=>36.6751, "bb_pure"=>87.2728, "f_bb"=>0.420236, "f_sc"=>0.688274),
    "PRO" => OrderedDict("n"=>97983, "sc"=>130.021, "sc_pure"=>192.851, "bb"=>39.646, "bb_pure"=>91.4236, "f_bb"=>0.433651, "f_sc"=>0.674206),
    "MET" => OrderedDict("n"=>36719, "sc"=>161.456, "sc_pure"=>222.969, "bb"=>39.909, "bb_pure"=>87.8564, "f_bb"=>0.454253, "f_sc"=>0.724115),
    "TRP" => OrderedDict("n"=>31497, "sc"=>230.509, "sc_pure"=>291.634, "bb"=>37.4479, "bb_pure"=>88.0122, "f_bb"=>0.425485, "f_sc"=>0.790406),
    "GLY" => OrderedDict("n"=>151702, "sc"=>0.0, "sc_pure"=>0.0, "bb"=>87.4824, "bb_pure"=>87.4824, "f_bb"=>1.0, "f_sc"=>1.0),
    "SER" => OrderedDict("n"=>129971, "sc"=>85.2744, "sc_pure"=>148.711, "bb"=>45.9041, "bb_pure"=>87.7078, "f_bb"=>0.523376, "f_sc"=>0.573423),
    "THR" => OrderedDict("n"=>118719, "sc"=>117.741, "sc_pure"=>179.107, "bb"=>39.4432, "bb_pure"=>87.0271, "f_bb"=>0.453229, "f_sc"=>0.657378),
    "TYR" => OrderedDict("n"=>78202, "sc"=>202.043, "sc_pure"=>263.305, "bb"=>38.86, "bb_pure"=>87.3211, "f_bb"=>0.445025, "f_sc"=>0.767335),
    "GLN" => OrderedDict("n"=>82612, "sc"=>156.603, "sc_pure"=>218.196, "bb"=>40.2219, "bb_pure"=>88.4311, "f_bb"=>0.454839, "f_sc"=>0.717718),
    "ASN" => OrderedDict("n"=>93411, "sc"=>129.155, "sc_pure"=>191.405, "bb"=>40.9267, "bb_pure"=>88.7074, "f_bb"=>0.461367, "f_sc"=>0.674776),
    "ASP" => OrderedDict("n"=>125837, "sc"=>121.985, "sc_pure"=>184.224, "bb"=>41.1593, "bb_pure"=>89.0758, "f_bb"=>0.46207, "f_sc"=>0.662158),
    "GLU" => OrderedDict("n"=>142356, "sc"=>149.532, "sc_pure"=>211.134, "bb"=>40.8443, "bb_pure"=>89.038, "f_bb"=>0.458729, "f_sc"=>0.708232),
    "HIS" => OrderedDict("n"=>51092, "sc"=>165.044, "sc_pure"=>226.829, "bb"=>40.172, "bb_pure"=>87.7198, "f_bb"=>0.457958, "f_sc"=>0.727611),
    "LYS" => OrderedDict("n"=>117537, "sc"=>178.395, "sc_pure"=>240.241, "bb"=>41.6446, "bb_pure"=>88.6705, "f_bb"=>0.469655, "f_sc"=>0.742567),
    "ARG" => OrderedDict("n"=>110038, "sc"=>210.225, "sc_pure"=>272.22, "bb"=>41.319, "bb_pure"=>88.3584, "f_bb"=>0.467629, "f_sc"=>0.772262),
    "CYS" => OrderedDict("n"=>29556, "sc"=>92.4815, "sc_pure"=>155.731, "bb"=>44.936, "bb_pure"=>87.0777, "f_bb"=>0.516044, "f_sc"=>0.593854),
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
    cosolvent = lowercase(cosolvent)
    col = cosolvent_column_Accessibility[cosolvent]
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