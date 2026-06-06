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
# These are the GTFEapp (apparent transfer free energies of backbones and side-chains - without corrections for
# glycine activity or specific side-chain activity)
#
const tfe_sc_bb_AutonBolenApp = Dict{String,NTuple{9,Float32}}(
#                TMAO   Sarcosine     Betaine     Proline    Sorbitol    Sucrose     UreaAPP    Glycerol  Trehalose
    "ALA" => ( -14.64,      10.91,       4.77,      -0.07,      16.57,      22.05,      -4.69,      7.76,     33.25),
    "PHE" => (  -9.32,     -12.64,    -112.93,     -71.26,      26.38,     -96.35,     -83.11,     59.77,    -17.88),
    "LEU" => (  11.62,      38.33,     -17.73,       4.77,      39.07,      37.11,     -54.57,    -34.42,     96.18),
    "ILE" => ( -25.43,      39.98,      -1.27,      -2.72,      36.90,      28.12,     -38.43,     36.23,     79.66),
    "VAL" => (  -1.02,      29.32,     -19.63,       7.96,      24.65,      33.92,     -21.65,     -1.37,     96.79),
    "PRO" => (-137.73,     -34.23,    -125.16,     -63.96,      -4.48,     -73.02,     -17.65,    -60.55,    -94.67),
    "MET" => (  -7.65,       8.18,     -14.16,     -35.12,      20.97,      -6.66,     -48.34,     13.87,     29.19),
    "TRP" => (-152.87,    -113.03,    -369.93,    -198.37,     -67.23,    -215.27,    -141.46,   -122.65,   -206.30),
    "GLY" => (      0,          0,          0,          0,          0,          0,       0.00,      0.00,      0.00),
    "SER" => ( -39.04,     -27.98,     -41.85,     -33.49,      -1.58,      -2.79,     -20.56,      6.31,     -0.98),
    "THR" => (   3.57,      -7.54,       0.33,     -18.33,      13.20,      20.82,     -22.09,     17.53,     26.32),
    "TYR" => (-114.32,     -26.37,    -213.09,    -138.41,     -53.50,     -78.41,     -45.08,   -149.50,    -80.32),
    "GLN" => (  41.41,     -10.19,       7.57,     -32.26,     -23.98,     -40.87,     -54.81,     -2.76,    -36.34),
    "ASN" => (  55.69,     -40.93,      33.17,     -17.71,     -21.21,     -28.28,     -38.79,     51.57,     48.67),
    "ASP" => ( -66.67,     -14.20,    -116.56,     -90.51,     -83.88,     -37.17,       3.55,    -85.46,    -96.54),
    "GLU" => ( -83.25,     -12.61,    -112.08,     -89.17,     -70.05,     -41.65,       0.62,    -74.20,    -85.92),
    "HIS" => (  42.07,     -20.80,     -35.97,     -45.10,     -42.45,    -118.66,     -50.51,    -17.17,    -98.75),
    "LYS" => (-110.23,     -27.42,    -171.99,     -59.87,     -32.47,     -39.60,     -22.76,    -34.01,    -50.07),
    "ARG" => (-109.27,     -32.24,    -109.45,     -60.18,     -24.65,     -79.32,     -21.17,    -30.74,    -50.33),
    "CYS" => (      0,          0,          0,          0,          0,          0,       0.00,      0.00,      0.00), # not reported
    "BB"  => (     90,         52,         67,         48,         35,         62,        -39,        14,        62),
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
# Glycine activity corrections
#
const γG = Dict{String,Float32}(
    "water" => 0.0,
    "tmao" => 0.0,
    "sarcosine" => 0.0,
    "betaine" => 0.0,
    "proline" => 0.0,
    "glycerol" => 0.0, # ~ -4
    "sorbitol" => 0.0, # ~ -0-2
    "sucrose" => 0.0,
    "trehalose" => 0.0,
    "urea" => -14.47,
)

const tfe_sc_bb_Accessibility = OrderedDict{String, NTuple{9, Float32}}(
  "ALA" => (40.5856, 61.0755, 60.652, 36.6777, 59.0326, 90.3858, 18.9989, 25.8121, 112.143),
  "PHE" => (48.1216, 18.2637, -103.547, -61.5004, 58.1674, -85.0998, -90.2562, 87.9847, 18.082),
  "LEU" => (80.1472, 90.1394, 23.0918, 40.7663, 79.0694, 95.5633, -55.5939, -37.7579, 177.457),
  "ILE" => (31.7858, 94.3353, 48.2037, 32.0187, 77.3769, 85.3133, -33.2934, 60.7842, 156.93),
  "VAL" => (67.3601, 82.9165, 22.3867, 48.435, 63.0398, 97.3343, -10.5619, 8.69635, 189.817),
  "PRO" => (-137.019, -10.4523, -136.29, -58.6876, 21.0876, -61.2657, -4.81861, -80.6002, -94.0717),
  "MET" => (49.9731, 46.0709, 25.6485, -15.7137, 52.1161, 32.5287, -46.2286, 28.3287, 81.4598),
  "TRP" => (-132.359, -107.508, -420.274, -217.362, -61.2466, -229.309, -159.522, -144.787, -218.041),
  "GLY" => (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.47, 0.0, 0.0),
  "SER" => (0.667515, -8.83995, -21.254, -21.1938, 23.3974, 41.4619, -10.3593, 21.1675, 44.5407),
  "THR" => (74.0381, 28.2891, 51.6053, 8.95793, 46.6141, 78.7066, -11.4971, 37.1284, 87.005),
  "TYR" => (-88.3386, 0.0347134, -230.495, -146.787, -45.9923, -60.413, -39.3987, -183.277, -62.8714),
  "GLN" => (117.767, 21.0795, 55.655, -12.0026, -9.38114, -14.4033, -55.5322, 5.65249, -8.16726),
  "ASN" => (145.884, -22.1609, 96.5345, 8.76978, -5.68575, 3.45808, -35.4121, 85.1697, 115.504),
  "ASP" => (-31.9228, 17.1044, -121.592, -97.3485, -97.1226, -9.18955, 26.3394, -114.722, -95.9693),
  "GLU" => (-54.2385, 17.6901, -109.595, -90.7751, -73.1392, -15.6754, 20.8462, -93.0514, -76.8325),
  "HIS" => (117.87, 6.57968, -4.02962, -29.2524, -34.3753, -120.174, -49.1562, -14.0093, -93.018),
  "LYS" => (-91.1235, -3.88577, -188.748, -50.0522, -21.4577, -13.9195, -11.146, -36.8474, -27.9965),
  "ARG" => (-85.7403, -9.7253, -100.048, -48.1802, -10.3476, -64.2786, -8.63445, -31.0476, -26.9185),
  "CYS" => (64.8481, 37.4678, 48.2758, 34.5856, 25.2187, 44.6731, 22.5536, 10.0875, 44.6731),
  "BB"  => (90.0, 52.0, 67.0, 48.0, 35.0, 62.0, -39.0, 14.0, 62.0),
)

#
# Backbone accessibility parameter - [-Inf,1]
#
# if acc ==  1 -> the accessibility is 1.0 (full access)
# if acc ==  0 -> the accessibility is defined by the ratio of ASAs of actual and pure BBs
# if acc == -1 -> the accessibility is zero (full shielding)
#
_acc_f(x, a) = x + (a + abs(a) * (1-2x))/2
#
_acc(x) = Dict{String,Float32}(
  "ALA" => x, "PHE" => x, "LEU" => x, "ILE" => x, "VAL" => x,
  "PRO" => x, "MET" => x, "TRP" => x, "GLY" => x, "SER" => x,
  "THR" => x, "TYR" => x, "GLN" => x, "ASN" => x, "ASP" => x,
  "GLU" => x, "HIS" => x, "LYS" => x, "ARG" => x, "CYS" => x,
)

const acc_bb = Dict{String,Dict{String,Float32}}(
    "tmao" => _acc(0.f0),
    "sarcosine" => _acc(0.f0),
    "betaine" => _acc(0.f0),
    "proline" => _acc(0.f0),
    "glycerol" => _acc(0.f0),
    "sorbitol" => _acc(0.f0),
    "sucrose" => _acc(0.f0),
    "trehalose" => _acc(0.f0),
    "urea" => _acc(1.f0),
)

const acc_sc = Dict{String,Dict{String,Float32}}(
    "tmao" => _acc(0.f0),
    "sarcosine" => _acc(0.f0),
    "betaine" => _acc(0.f0),
    "proline" => _acc(0.f0),
    "glycerol" => _acc(0.f0),
    "sorbitol" => _acc(0.f0),
    "sucrose" => _acc(0.f0),
    "trehalose" => _acc(0.f0),
    "urea" => _acc(0.f0),
)

function set_acc()
    acc_bb["tmao"]["THR"] = 0.0 
    acc_bb["tmao"]["PHE"] = 0.0
    acc_bb["tmao"]["TRP"] = 0.0
    acc_bb["tmao"]["TYR"] = 0.0
end

function model_combination_rule(::Type{Accessibility}, cosolvent, restype)
    # voltar
    #tfe_sc_bb = tfe_sc_bb_Accessibility
    tfe_sc_bb = _tfe_sc_bb_Accessibility()
    col = cosolvent_column_Accessibility[lowercase(cosolvent)]
    # united model: all bb ASA contributions are the same
    bb_contribution = tfe_sc_bb["BB"][col] / f_acc["GLY"]["bb_pure"]
    sc_contribution = if restype == "GLY"
        0.0f0
    else
        tfe_sc_bb[restype][col] / f_acc[restype]["sc_pure"]
    end
    return bb_contribution, sc_contribution
end

#=

unused data / development

=#

#
# Transfer free energy of glycine in each cosolvent
#
# In molarity scale - cal/mol
# Data computed from the solubility table of: https://pubs.acs.org/doi/10.1021/bi035908r
#
const tfe_gly = Dict{String,Float32}(
    "water" => 0.0,
    "tmao" => 182.28,
    "sarcosine" => 59.84,
    "betaine" => 167.54,
    "proline" => 102.80,
    "glycerol" => 56.18,
    "sorbitol" => 54.08,
    "sucrose" => 144.48,
    "trehalose" => 143.50,
    "urea" => 18.74,
)

#
# This function generates the tfe_sc_bb_Accessibility dictionary from the data
#
function _tfe_sc_bb_Accessibility()
    set_acc()
    cs = collect(keys(cosolvent_column_Accessibility))
    data = OrderedDict{String,NTuple{9,Float32}}(
        [ 
            aa => ntuple(
                i -> begin 
                    f_sc = _acc_f(f_acc[aa]["f_sc"], acc_sc[cs[i]][aa])
                    f_bb = _acc_f(f_acc[aa]["f_bb"], acc_bb[cs[i]][aa])
                    GTFEapp = tfe_sc_bb_AutonBolenApp[aa][i] # apparent free energy
                    tfe_bb = tfe_sc_bb_AutonBolenApp["BB"][i]
                    γG_val = γG[cs[i]] # -14.47 for urea
                    tfe_sc = inv(f_sc) * (GTFEapp - γG_val + tfe_bb * (1 - f_bb))
                end,
            9)
            for aa in keys(f_acc)  
        ]...
    )
    data["BB"] = tfe_sc_bb_AutonBolenApp["BB"]
    return data
end
