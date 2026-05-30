#=


=#

export UniversalPure
struct UniversalPure <: MValueModel end

# Do not user underscores (_) in the following names:
const cosolvent_column_UniversalPure = OrderedDict(
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
# n: number of entries in the database
#
bb_sc_exp() = 
Dict{String, Tuple{Int64, Vararg{Float32, 6}}}(
#            1      2        3       4       5         6          7
#            n     bb     bb_pure    sc     sc_pure   f_bb       f_sc
  "GLN" => (240, 44.2294, 85.8413, 151.59, 208.669, 0.515246, 0.726461),
  "LYS" => (401, 45.4951, 86.9512, 166.6, 224.09, 0.523226, 0.743451),
  "GLY" => (583, 86.5575, 86.5575, 0.0, 0.0, 1.0, 1.0),
  "ASN" => (318, 43.7799, 87.2324, 127.355, 184.834, 0.501877, 0.689026),
  "TRP" => (100, 40.6636, 87.1785, 222.516, 279.354, 0.466441, 0.796537),
  "THR" => (400, 42.5826, 85.522, 112.726, 169.425, 0.497914, 0.665346),
  "VAL" => (459, 41.4824, 85.96, 121.062, 177.501, 0.482578, 0.682033),
  "SER" => (502, 48.8819, 87.0889, 83.2252, 141.471, 0.561288, 0.588283),
  "HIS" => (155, 42.671, 85.0285, 159.036, 216.355, 0.501843, 0.735071),
  "PRO" => (281, 42.7962, 89.859, 116.505, 176.266, 0.47626, 0.660964),
  "PHE" => (239, 41.5454, 85.956, 181.218, 238.035, 0.483333, 0.761309),
  "ASP" => (399, 44.0283, 87.96, 124.616, 182.095, 0.500549, 0.684344),
  "ILE" => (338, 40.4473, 86.4206, 146.245, 202.6, 0.468029, 0.721841),
  "TYR" => (299, 42.7139, 85.8969, 198.902, 255.919, 0.497269, 0.777209),
  "ARG" => (230, 45.136, 86.3419, 200.409, 257.924, 0.522759, 0.777008),
  "LEU" => (485, 42.5356, 86.8573, 145.019, 200.877, 0.489719, 0.721929),
  "ALA" => (510, 53.0357, 87.6011, 62.075, 120.762, 0.605423, 0.514025),
  "MET" => (116, 44.8261, 86.0533, 155.812, 212.591, 0.520911, 0.73292),
  "CYS" => (156, 46.4347, 86.3576, 103.173, 160.629, 0.537703, 0.642311),
  "GLU" => (341, 45.2395, 87.0802, 148.887, 206.045, 0.519516, 0.722597),
)

#
# Glycine activity corrections
#
γG() = Dict{String,Float32}(
    "water" => 0.0,
    "tmao" => 0.0,
    "sarcosine" => 0.0,
    "betaine" => 0.0,
    "proline" => 0.0,
    "glycerol" => 0.0,
    "sorbitol" => 0.0,
    "sucrose" => 0.0,
    "trehalose" => 0.0,
    "urea" => 14.47,
)

#
# Backbone accessibility parameter
#
acc() = Dict{String,Float32}(
    "tmao" => 1.0,
    "sarcosine" => 1.0,
    "betaine" => 1.0,
    "proline" => 1.0,
    "glycerol" => 1.0,
    "sorbitol" => 1.0,
    "sucrose" => 1.0,
    "trehalose" => 1.0,
    "urea" => 0.0,
)

function tfe_sc_bb_UniversalPure()
    cs = collect(keys(cosolvent_column_UniversalPure))
    data = Dict{String,NTuple{9,Float32}}(
        [ 
            aa => ntuple(
                i -> begin 
                    f_sc = bb_sc_exp()[aa][7]
                    f_bb = bb_sc_exp()[aa][6]
                    _acc = 1 - f_bb * acc()[cs[i]]
                    GTFEapp = tfe_sc_bb_AutonBolenApp[aa][i] # apparent free energy
                    f_bb_gly = bb_sc_exp()["GLY"][6]  # = 1
                    tfe_bb = tfe_sc_bb_AutonBolenApp["BB"][i]
                    γG_val = γG()[cs[i]]              # e.g., +14.47 for urea
                    tfe_sc = (1/f_sc) * (GTFEapp + γG_val + tfe_bb * (f_bb_gly - f_bb * _acc ))
                end,
            9)
            for aa in keys(bb_sc_exp())  
        ]...
    )
    data["BB"] = ( 90, 52, 67, 48, 35, 62, -39, 14, 62)
    return data
end

function model_combination_rule(::Type{UniversalPure}, cosolvent, restype)
    tfe_sc_bb = tfe_sc_bb_UniversalPure()
    col = cosolvent_column_UniversalPure[lowercase(cosolvent)]
    # united model: all bb ASA contributions are the same
    bb_contribution = tfe_sc_bb["BB"][col] / bb_sc_exp()["GLY"][3]
    sc_contribution = if restype == "GLY"
        0.0f0
    else
        tfe_sc_bb[restype][col] / bb_sc_exp()[restype][5]
    end
    return bb_contribution, sc_contribution
end

#=

unused data

=#

#
# Transfer free energy of glycine in each cosolvent
#
# In molarity scale - cal/mol
# Data computed from the solubility table of: https://pubs.acs.org/doi/10.1021/bi035908r
#
tfe_gly() = Dict{String,Float32}(
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
