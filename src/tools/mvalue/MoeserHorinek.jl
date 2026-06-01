#=

Amino acid side-chain and peptide backbone unit transfer free energies (cal/mol) from water to 1M osmolyte
Values for Urea GTFE+ values from Table S1 of https://doi.org/10.1021/jp409934q (https://pubs.acs.org/doi/10.1021/jp409934q#_i14)

=#
export MoeserHorinek
struct MoeserHorinek <: MValueModel end

# Do not user underscores (_) in the following names:
const cosolvent_column_MoeserHorinek = Dict(
    "tmao" => 1,
    "sarcosine"=> 2,
    "betaine"=> 3,
    "proline" => 4,
    "sorbitol" => 5,
    "sucrose" => 6,
    "urea" => 7,
    "glycerol" => 8,
    "trehalose" => 9,
)

#
# Note: Urea is the the GTFE+ pf the MH paper. 
#       For all other osmolytes, the values are GTFEapp, because there's not reliable Gly-activity 
#       correction (which, at the same time, appear to be small in any case for the data existing).
#
const tfe_sc_bb_MoeserHorinek = Dict{String,NTuple{9,Float32}}(
#                TMAO   Sarcosine     Betaine     Proline    Sorbitol    Sucrose    Urea    Glycerol  Trehalose
    "ALA" => ( -14.64,      10.91,       4.77,      -0.07,      16.57,      22.05,    1.01,     7.76,     33.25),
    "PHE" => (  -9.32,     -12.64,    -112.93,     -71.26,      26.38,     -96.35,  -68.64,    59.77,    -17.88),
    "LEU" => (  11.62,      38.33,     -17.73,       4.77,      39.07,      37.11,  -40.10,   -34.42,     96.18),
    "ILE" => ( -25.43,      39.98,      -1.27,      -2.72,      36.90,      28.12,  -23.96,    36.23,     79.66),
    "VAL" => (  -1.02,      29.32,     -19.63,       7.96,      24.65,      33.92,   -7.18,    -1.37,     96.79),
    "PRO" => (-137.73,     -34.23,    -125.16,     -63.96,      -4.48,     -73.02,   -3.18,   -60.55,    -94.67),
    "MET" => (  -7.65,       8.18,     -14.16,     -35.12,      20.97,      -6.66,  -33.87,    13.87,     29.19),
    "TRP" => (-152.87,    -113.03,    -369.93,    -198.37,     -67.23,    -215.27, -126.90,  -122.65,   -206.30),
    "GLY" => (      0,          0,          0,          0,          0,          0,    0.00,     0.00,      0.00),
    "SER" => ( -39.04,     -27.98,     -41.85,     -33.49,      -1.58,      -2.79,   -6.09,     6.31,     -0.98),
    "THR" => (   3.57,      -7.54,       0.33,     -18.33,      13.20,      20.82,   -7.62,    17.53,     26.32),
    "TYR" => (-114.32,     -26.37,    -213.09,    -138.41,     -53.50,     -78.41,  -30.61,  -149.50,    -80.32),
    "GLN" => (  41.41,     -10.19,       7.57,     -32.26,     -23.98,     -40.87,  -40.34,    -2.76,    -36.34),
    "ASN" => (  55.69,     -40.93,      33.17,     -17.71,     -21.21,     -28.28,  -24.32,    51.57,     48.67),
    "ASP" => ( -66.67,     -14.20,    -116.56,     -90.51,     -83.88,     -37.17,   18.02,   -85.46,    -96.54),
    "GLU" => ( -83.25,     -12.61,    -112.08,     -89.17,     -70.05,     -41.65,   15.09,   -74.20,    -85.92),
    "HIS" => (  42.07,     -20.80,     -35.97,     -45.10,     -42.45,    -118.66,  -36.04,   -17.17,    -98.75),
    "LYS" => (-110.23,     -27.42,    -171.99,     -59.87,     -32.47,     -39.60,   -8.29,   -34.01,    -50.07),
    "ARG" => (-109.27,     -32.24,    -109.45,     -60.18,     -24.65,     -79.32,   -6.70,   -30.74,    -50.33),
    "CYS" => (      0,          0,          0,          0,          0,          0,    0.00,     0.00,      0.00), # not reported
    "BB"  => (     90,         52,         67,         48,         35,         62,     -39,       14,        62),
)

# united model: all bb ASA contributions are the same
function model_combination_rule(::Type{MoeserHorinek}, cosolvent, restype)
    tfe_sc_bb = tfe_sc_bb_MoeserHorinek
    col = cosolvent_column_MoeserHorinek[lowercase(cosolvent)]
    bb_contribution = tfe_sc_bb["BB"][col] / first(isolated_ASA["GLY"])
    sc_contribution = if restype == "GLY"
        0.0f0
    else
        tfe_sc_bb[restype][col] / last(isolated_ASA[restype])
    end
    return bb_contribution, sc_contribution
end