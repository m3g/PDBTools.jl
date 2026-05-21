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
)
cosolvent_column(::Type{AutonBolen}) = cosolvent_column_AutonBolen

const tfe_sc_bb_AutonBolen = Dict{String,NTuple{8,Float32}}(
#                TMAO   Sarcosine     Betaine     Proline    Sorbitol    Sucrose      UreaWrong   UreaAPP 
    "ALA" => ( -14.64,      10.91,       4.77,      -0.07,      16.57,      22.05,       0.63,       -4.69),
    "PHE" => (  -9.32,     -12.64,    -112.93,     -71.26,      26.38,     -96.35,     -42.84,      -83.11),
    "LEU" => (  11.62,      38.33,     -17.73,       4.77,      39.07,      37.11,     -14.30,      -54.57),
    "ILE" => ( -25.43,      39.98,      -1.27,      -2.72,      36.90,      28.12,       1.84,      -38.43),
    "VAL" => (  -1.02,      29.32,     -19.63,       7.96,      24.65,      33.92,      18.62,      -21.65),
    "PRO" => (-137.73,     -34.23,    -125.16,     -63.96,      -4.48,     -73.02,      22.62,      -17.65),
    "MET" => (  -7.65,       8.18,     -14.16,     -35.12,      20.97,      -6.66,      -8.07,      -48.34),
    "TRP" => (-152.87,    -113.03,    -369.93,    -198.37,     -67.23,    -215.27,    -101.19,     -141.46),
    "GLY" => (      0,          0,          0,          0,          0,          0,       0.00,           0),
    "SER" => ( -39.04,     -27.98,     -41.85,     -33.49,      -1.58,      -2.79,      19.71,      -20.56),
    "THR" => (   3.57,      -7.54,       0.33,     -18.33,      13.20,      20.82,      18.18,      -22.09),
    "TYR" => (-114.32,     -26.37,    -213.09,    -138.41,     -53.50,     -78.41,      -4.81,      -45.08),
    "GLN" => (  41.41,     -10.19,       7.57,     -32.26,     -23.98,     -40.87,     -14.54,      -54.81),
    "ASN" => (  55.69,     -40.93,      33.17,     -17.71,     -21.21,     -28.28,       1.48,      -38.79),
    "ASP" => ( -66.67,     -14.20,    -116.56,     -90.51,     -83.88,     -37.17,      43.82,        3.55),
    "GLU" => ( -83.25,     -12.61,    -112.08,     -89.17,     -70.05,     -41.65,      40.89,        0.62),
    "HIS" => (  42.07,     -20.80,     -35.97,     -45.10,     -42.45,    -118.66,     -10.24,      -50.51),
    "LYS" => (-110.23,     -27.42,    -171.99,     -59.87,     -32.47,     -39.60,      17.51,      -22.76),
    "ARG" => (-109.27,     -32.24,    -109.45,     -60.18,     -24.65,     -79.32,      19.10,      -21.17),
    "CYS" => (      0,          0,          0,          0,          0,          0,       0.00,           0), # not reported
    "BB"  => (     90,         52,         67,         48,         35,         62,        -39,         -39),
)
tfe_sc_bb(::Type{AutonBolen}) = tfe_sc_bb_AutonBolen

#
# MoserHorinekFit consists of parameters to be used in the Moeser&Horinek model, but with Gly transfer free energy
# corrections that are obtained by fitting the predicted m-values of the Auton&Bolen model.
#

# Do not user underscores (_) in the following names:
const cosolvent_column_MoeserHorinekFit = Dict(
    "tmao" => 1,
    "sarcosine"=> 2,    
    "betaine"=> 3,    
    "proline" => 4,    
    "sorbitol" => 5,   
    "sucrose" => 6,
    "urea" => 7,
)
cosolvent_column(::Type{MoeserHorinekFit}) = cosolvent_column_MoeserHorinekFit

# These where obtained by minmizing the sum of squared residues relative to AutonBolen predictions
# using: BlackBoxOptim.jl - BlackBoxOptim.AdaptiveDiffEvoRandBin{3}
#               TMAO  Sarcosine   Betaine    Proline    Sorbitol    Sucrose       Urea
const γG =   ( 47.74,     27.57,    35.57,     25.44,      18.57,     32.91,     14.74)

const tfe_sc_bb_MoeserHorinekFit = Dict{String, NTuple{7,Float32}}(
#                TMAO    Sarcosine     Betaine     Proline    Sorbitol    Sucrose       Urea 
    "ALA" => ( -14.64,       10.91,       4.77,      -0.07,      16.57,      22.05,     -4.69) .+ γG,
    "PHE" => (  -9.32,      -12.64,    -112.93,     -71.26,      26.38,     -96.35,    -83.11) .+ γG,
    "LEU" => (  11.62,       38.33,     -17.73,       4.77,      39.07,      37.11,    -54.57) .+ γG,
    "ILE" => ( -25.43,       39.98,      -1.27,      -2.72,      36.90,      28.12,    -38.43) .+ γG,
    "VAL" => (  -1.02,       29.32,     -19.63,       7.96,      24.65,      33.92,    -21.65) .+ γG,
    "PRO" => (-137.73,      -34.23,    -125.16,     -63.96,      -4.48,     -73.02,    -17.65) .+ γG,
    "MET" => (  -7.65,        8.18,     -14.16,     -35.12,      20.97,      -6.66,    -48.34) .+ γG,
    "TRP" => (-152.87,     -113.03,    -369.93,    -198.37,     -67.23,    -215.27,   -141.46) .+ γG,
    "GLY" => (      0,           0,          0,          0,          0,          0,         0),
    "SER" => ( -39.04,      -27.98,     -41.85,     -33.49,      -1.58,      -2.79,    -20.56) .+ γG,
    "THR" => (   3.57,       -7.54,       0.33,     -18.33,      13.20,      20.82,    -22.09) .+ γG,
    "TYR" => (-114.32,      -26.37,    -213.09,    -138.41,     -53.50,     -78.41,    -45.08) .+ γG,
    "GLN" => (  41.41,      -10.19,       7.57,     -32.26,     -23.98,     -40.87,    -54.81) .+ γG,
    "ASN" => (  55.69,      -40.93,      33.17,     -17.71,     -21.21,     -28.28,    -38.79) .+ γG,
    "ASP" => ( -66.67,      -14.20,    -116.56,     -90.51,     -83.88,     -37.17,      3.55) .+ γG,
    "GLU" => ( -83.25,      -12.61,    -112.08,     -89.17,     -70.05,     -41.65,      0.62) .+ γG,
    "HIS" => (  42.07,      -20.80,     -35.97,     -45.10,     -42.45,    -118.66,    -50.51) .+ γG,
    "LYS" => (-110.23,      -27.42,    -171.99,     -59.87,     -32.47,     -39.60,    -22.76) .+ γG,
    "ARG" => (-109.27,      -32.24,    -109.45,     -60.18,     -24.65,     -79.32,    -21.17) .+ γG,
    "CYS" => (   0.00,           0,          0,          0,          0,          0,         0),
    "BB"  => (     90,          52,         67,         48,         35,         62,       -39),
)
tfe_sc_bb(::Type{MoeserHorinekFit}) = tfe_sc_bb_MoeserHorinekFit

