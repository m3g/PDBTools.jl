#
# Data section
#

#= 

Isolated ASA values are from the Supporting Table 2 of https://doi.org/10.1073/pnas.0507053102
(https://www.pnas.org/doi/suppl/10.1073/pnas.0507053102/suppl_file/07053table2.pdf)

=#
const isolated_ASA = Dict{String,Tuple{Float32,Float32}}( 
                # BB      SC   (Å^2)
    "ALA"	=> (46.2,	71.9),
    "PHE"	=> (38.4,	184.4),
    "LEU"	=> (35.3,	157.8),
    "ILE"	=> (30.9,	150.1),
    "VAL"	=> (36.1,	128.4),
    "PRO"	=> (35.6,	111.0),
    "MET"	=> (38.6,	164.8),
    "TRP"	=> (37.4,	228.9),
    "GLY"	=> (88.1,	0),
    "SER"	=> (44.0,	85.8),
    "THR"	=> (37.9,	114.6),
    "TYR"	=> (38.7,	198.1),
    "GLN"	=> (37.8,	155.4),
    "ASN"	=> (40.2,	125.3),
    "ASP"	=> (40.5,	118.2),
    "GLU"	=> (37.8,	148.4),
    "HIS"	=> (40.4,	162.1),
    "LYS"	=> (38.7,	187.1),
    "ARG"	=> (39.1,	216.9),
    "CYS"	=> (42.6,	103.5),
) 

#=

Amino acid side-chain and peptide backbone unit transfer free energies (cal/mol) from water to 1M osmolyte
Values for Urea GTFE+ values from Table S1 of https://doi.org/10.1021/jp409934q (https://pubs.acs.org/doi/10.1021/jp409934q#_i14)

Only Urea values are available for this model for now.

=#

# Do not user underscores (_) in the following names:
const cosolvent_column_MoeserHorinek = Dict(
    "urea" => 1,
)
cosolvent_column(::Type{MoeserHorinek}) = cosolvent_column_MoeserHorinek

const tfe_sc_bb_MoeserHorinek = Dict{String, NTuple{1,Float32}}(
#               Urea 
    "ALA" => (   1.01, ),
    "PHE" => ( -68.64, ),
    "LEU" => ( -40.10, ),
    "ILE" => ( -23.96, ),
    "VAL" => (  -7.18, ),
    "PRO" => (  -3.18, ),
    "MET" => ( -33.87, ),
    "TRP" => (-126.90, ),
    "GLY" => (   0.00, ),
    "SER" => (  -6.09, ),
    "THR" => (  -7.62, ),
    "TYR" => ( -30.61, ),
    "GLN" => ( -40.34, ),
    "ASN" => ( -24.32, ),
    "ASP" => (  18.02, ),
    "GLU" => (  15.09, ),
    "HIS" => ( -36.04, ),
    "LYS" => (  -8.29, ),
    "ARG" => (  -6.70, ),
    "CYS" => (   0.00, ),
    "BB"  => (    -39, ),
)
tfe_sc_bb(::Type{MoeserHorinek}) = tfe_sc_bb_MoeserHorinek

#=

Amino acid side-chain and peptide backbone unit transfer free energies (cal/mol) from water to 1M osmolyte
Supplementary Table 1 of https://doi.org/10.1073/pnas.0507053102

UreaWrong from GTFE* from Supplementary Table S1 of Moeser and Horinek and originally from https://doi.org/10.1073/pnas.0706251104

The "urea" column selection points to "UreaWrong" which is what is output from the server.

=#

# Do not user underscores (_) in the following names:
const cosolvent_column_AutonBolen = Dict(
    "tmao" => 1,
    "sarcosine"=> 2,
    "betaine"=> 3,
    "proline" => 4,
    "sorbitol" => 5,
    "sucrose" => 6,
    "urea" => 7,
    "urea-app" => 8,
    "glycerol" => 9,
    "trehalose" => 10,
)
cosolvent_column(::Type{AutonBolen}) = cosolvent_column_AutonBolen

const tfe_sc_bb_AutonBolen = Dict{String,NTuple{10,Float32}}(
#                TMAO   Sarcosine     Betaine     Proline    Sorbitol    Sucrose      UreaWrong   UreaAPP    Glycerol  Trehalose
    "ALA" => ( -14.64,      10.91,       4.77,      -0.07,      16.57,      22.05,       0.63,       -4.69,      7.76,     33.25),
    "PHE" => (  -9.32,     -12.64,    -112.93,     -71.26,      26.38,     -96.35,     -42.84,      -83.11,     59.77,    -17.88),
    "LEU" => (  11.62,      38.33,     -17.73,       4.77,      39.07,      37.11,     -14.30,      -54.57,    -34.42,     96.18),
    "ILE" => ( -25.43,      39.98,      -1.27,      -2.72,      36.90,      28.12,       1.84,      -38.43,     36.23,     79.66),
    "VAL" => (  -1.02,      29.32,     -19.63,       7.96,      24.65,      33.92,      18.62,      -21.65,     -1.37,     96.79),
    "PRO" => (-137.73,     -34.23,    -125.16,     -63.96,      -4.48,     -73.02,      22.62,      -17.65,    -60.55,    -94.67),
    "MET" => (  -7.65,       8.18,     -14.16,     -35.12,      20.97,      -6.66,      -8.07,      -48.34,     13.87,     29.19),
    "TRP" => (-152.87,    -113.03,    -369.93,    -198.37,     -67.23,    -215.27,    -101.19,     -141.46,   -122.65,   -206.30),
    "GLY" => (      0,          0,          0,          0,          0,          0,       0.00,        0.00,      0.00,      0.00),
    "SER" => ( -39.04,     -27.98,     -41.85,     -33.49,      -1.58,      -2.79,      19.71,      -20.56,      6.31,     -0.98),
    "THR" => (   3.57,      -7.54,       0.33,     -18.33,      13.20,      20.82,      18.18,      -22.09,     17.53,     26.32),
    "TYR" => (-114.32,     -26.37,    -213.09,    -138.41,     -53.50,     -78.41,      -4.81,      -45.08,   -149.50,    -80.32),
    "GLN" => (  41.41,     -10.19,       7.57,     -32.26,     -23.98,     -40.87,     -14.54,      -54.81,     -2.76,    -36.34),
    "ASN" => (  55.69,     -40.93,      33.17,     -17.71,     -21.21,     -28.28,       1.48,      -38.79,     51.57,     48.67),
    "ASP" => ( -66.67,     -14.20,    -116.56,     -90.51,     -83.88,     -37.17,      43.82,        3.55,    -85.46,    -96.54),
    "GLU" => ( -83.25,     -12.61,    -112.08,     -89.17,     -70.05,     -41.65,      40.89,        0.62,    -74.20,    -85.92),
    "HIS" => (  42.07,     -20.80,     -35.97,     -45.10,     -42.45,    -118.66,     -10.24,      -50.51,    -17.17,    -98.75),
    "LYS" => (-110.23,     -27.42,    -171.99,     -59.87,     -32.47,     -39.60,      17.51,      -22.76,    -34.01,    -50.07),
    "ARG" => (-109.27,     -32.24,    -109.45,     -60.18,     -24.65,     -79.32,      19.10,      -21.17,    -30.74,    -50.33),
    "CYS" => (      0,          0,          0,          0,          0,          0,       0.00,        0.00,      0.00,      0.00), # not reported
    "BB"  => (     90,         52,         67,         48,         35,         62,        -39,         -39,        14,        62),
)
tfe_sc_bb(::Type{AutonBolen}) = tfe_sc_bb_AutonBolen

#=

Here we have the apparent GTFEs reported by Auton and Bolen in
Supplementary Table 1 of https://doi.org/10.1073/pnas.0507053102

However, we apply the universal backbone description of Moeser and Horinek, but without the error 
compensation of Gly-activity that they applied. 

# Compute proper universal backbone side-chain TFEs
# GTFE^UB_aa = GTFE^est_aa + TFE_bb × (1 - ASA^bb_aa / 88.1)

=#

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
cosolvent_column(::Type{MoeserHorinekApp}) = cosolvent_column_MoeserHorinekApp

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
tfe_sc_bb(::Type{MoeserHorinekApp}) = tfe_sc_bb_MoeserHorinekApp