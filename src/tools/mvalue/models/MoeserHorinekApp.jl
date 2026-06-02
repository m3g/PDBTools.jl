#=

Here we have the apparent GTFEs reported by Auton and Bolen in
Supplementary Table 1 of https://doi.org/10.1073/pnas.0507053102

However, we apply the universal backbone description of Moeser and Horinek, but without the error 
compensation of Gly-activity that they applied. 

# Compute proper universal backbone side-chain TFEs
# GTFE^UB_aa = GTFE^est_aa + TFE_bb × (1 - ASA^bb_aa / 88.1)

=#

export MoeserHorinekApp
struct MoeserHorinekApp <: MValueModel end

# Do not user underscores (_) in the following names:
const cosolvent_column_MoeserHorinekApp = Dict(
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

const tfe_sc_bb_MoeserHorinekApp = Dict{String,NTuple{9,Float32}}(
#                TMAO   Sarcosine     Betaine     Proline    Sorbitol    Sucrose    UreaAPP    Glycerol  Trehalose
    "ALA" => (   28.16,      35.64,      36.63,      22.76,      33.22,      51.54,    -23.24,     14.42,     62.74),
    "PHE" => (   41.45,      16.69,     -75.13,     -44.18,      46.12,     -61.37,   -105.11,     67.67,     17.10),
    "LEU" => (   65.56,      69.49,      22.42,      33.54,      60.05,      74.27,    -77.94,    -26.03,    133.34),
    "ILE" => (   33.00,      73.74,      42.23,      28.44,      59.62,      68.37,    -63.75,     45.32,    119.91),
    "VAL" => (   52.10,      60.01,      19.92,      36.29,      45.31,      70.51,    -44.67,      6.89,    133.38),
    "PRO" => (  -84.10,      -3.24,     -85.23,     -35.36,      16.38,     -36.07,    -40.89,    -52.21,    -57.72),
    "MET" => (   42.92,      37.40,      23.48,      -8.15,      40.64,      28.18,    -70.25,     21.74,     64.03),
    "TRP" => ( -101.08,     -83.10,    -331.37,    -170.75,     -47.09,    -179.59,   -163.90,   -114.59,   -170.62),
    "GLY" => (    0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      0.00,      0.00,      0.00),
    "SER" => (    6.01,      -1.95,      -8.31,      -9.46,      15.94,      28.25,    -40.08,     13.32,     30.06),
    "THR" => (   54.85,      22.09,      38.51,       9.02,      33.14,      56.15,    -44.31,     25.51,     61.65),
    "TYR" => (  -63.85,       2.79,    -175.52,    -111.50,     -33.87,     -43.64,    -66.95,   -141.65,    -45.55),
    "GLN" => (   92.79,      19.50,      45.82,      -4.85,      -4.00,      -5.47,    -77.08,      5.23,     -0.94),
    "ASN" => (  104.62,     -12.66,      69.60,       8.39,      -2.18,       5.43,    -59.99,     59.18,     82.38),
    "ASP" => (  -18.04,      13.90,     -80.36,     -64.58,     -64.97,      -3.67,    -17.52,    -77.90,    -63.04),
    "GLU" => (  -31.87,      17.08,     -73.83,     -61.76,     -50.07,      -6.25,    -21.65,    -66.21,    -50.52),
    "HIS" => (   90.80,       7.35,       0.31,     -19.11,     -23.50,     -85.09,    -71.63,     -9.59,    -65.18),
    "LYS" => (  -59.76,       1.74,    -134.42,     -32.96,     -12.84,      -4.83,    -44.63,    -26.16,    -15.30),
    "ARG" => (  -59.21,      -3.32,     -72.19,     -33.48,      -5.18,     -44.84,    -42.86,    -22.95,    -15.85),
    "CYS" => (   46.48,      26.86,      34.60,      24.79,      18.08,      32.02,    -20.14,      7.23,     32.02), # not reported in App; correction applied to zero
    "BB"  => (      90,         52,         67,         48,         35,         62,       -39,        14,        62),
)

function model_combination_rule(::Type{MoeserHorinekApp}, cosolvent, restype)
    tfe_sc_bb = tfe_sc_bb_MoeserHorinekApp
    col = cosolvent_column_MoeserHorinekApp[lowercase(cosolvent)]
    bb_contribution = tfe_sc_bb["BB"][col] / first(isolated_ASA["GLY"])
    sc_contribution = if restype == "GLY"
        0.0f0
    else
        tfe_sc_bb[restype][col] / last(isolated_ASA[restype])
    end
    return bb_contribution, sc_contribution
end