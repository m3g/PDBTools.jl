#=

Amino acid side-chain and peptide backbone unit transfer free energies (cal/mol) from water to 1M osmolyte
Values for Urea GTFE+ values from Table S1 of https://doi.org/10.1021/jp409934q (https://pubs.acs.org/doi/10.1021/jp409934q#_i14)

Only Urea values are available for this model for now.

=#

export MoeserHorinek
struct MoeserHorinek <: MValueModel end
modelname(::Type{MoeserHorinek}) = "MoeserHorinek"
modeltype(::Type{MoeserHorinek}) = MoeserHorinek

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