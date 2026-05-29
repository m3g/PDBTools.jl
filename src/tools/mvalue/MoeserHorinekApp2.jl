#=


=#

export MoeserHorinekApp2
struct MoeserHorinekApp2 <: MValueModel end
modeltype(::Type{MoeserHorinekApp2}) = MoeserHorinek

# Do not user underscores (_) in the following names:
const cosolvent_column_MoeserHorinekApp2 = OrderedDict(
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
cosolvent_column(::Type{MoeserHorinekApp2}) = cosolvent_column_MoeserHorinekApp2

#
# Fractional bb exposure in the amino acid
# Compute by, for example:
#
# isolated_ASA["GLN"][1] / isolated_ASA["GLY"][1] 
#
# voltar: this is the crucial data to reconcile the models, we need better
# estimates of the protection that the side-chains promote to the backbones
#
bb_exp() = Dict{String,Float32}(
    "PHE" => 0.80, #0.435868,
    "GLN" => 0.80, #0.429058,
    "ASP" => 0.83, #0.459705,
    "LYS" => 0.73, #0.439274,
    "ILE" => 0.75, #0.350738,
    "TYR" => 0.74, #0.439274,
    "GLY" => 1.0, #1.0,
    "ASN" => 0.80, #0.4563,
    "ARG" => 0.77, #0.443814,
    "LEU" => 0.80, #0.400681,
    "TRP" => 0.74, #0.424518,
    "ALA" => 0.87, #0.524404,
    "THR" => 0.77, #0.430193,
    "VAL" => 0.77, #0.409762,
    "MET" => 0.7, #0.438138,
    "SER" => 0.79, #0.499432,
    "GLU" => 0.84, #0.429058,
    "PRO" => 0.74, #0.404086,
    "HIS" => 0.7, #0.45857,
    "CYS" => 0.7, #0.483541,
)

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

function tfe_sc_bb_MoeserHorinekApp2()
    cs = collect(keys(cosolvent_column_MoeserHorinekApp2))
    data = Dict{String,NTuple{9,Float32}}(
        [ 
            aa => ntuple(
                i -> begin 
                    tfe_sc_bb_AutonBolen[aa][i] + (1 - bb_exp()[aa]) * tfe_gly()[cs[i]] + bb_exp()[aa] * γG()[cs[i]] 
                end,
            9)
            for aa in keys(bb_exp())  
        ]...
    )
    data["BB"] = ( 90, 52, 67, 48, 35, 62, -39, 14, 62)
    return data
end

tfe_sc_bb(::Type{MoeserHorinekApp2}) = tfe_sc_bb_MoeserHorinekApp2()