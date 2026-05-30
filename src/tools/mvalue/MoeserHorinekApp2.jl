#=


=#

export MoeserHorinekApp2
struct MoeserHorinekApp2 <: MValueModel end
modeltype(::Type{MoeserHorinekApp2}) = MoeserHorinekApp2

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
bb_sc_exp() = 
Dict{String, Tuple{Int64, Vararg{Float32, 6}}}(
#            1      2        3       4       5         6          7
#            n     bb     bb_pure    sc     sc_pure   f_bb       f_sc
  "GLN" => (240, 44.2294, 85.8413, 144.595, 208.669, 0.515246, 0.692937),
  "LYS" => (401, 45.4951, 86.9512, 160.833, 224.09, 0.523226, 0.717716),
  "GLY" => (583, 86.5575, 86.5575,   0.0,     0.0, 1.0, 1.0),
  "ASN" => (318, 43.7799, 87.2324, 120.06,  184.834, 0.501877, 0.64956),
  "TRP" => (100, 40.6636, 87.1785, 210.652, 279.354, 0.466441, 0.754068),
  "THR" => (400, 42.5826, 85.522, 103.975,  169.425, 0.497914, 0.613693),
  "VAL" => (459, 41.4824, 85.96,  111.524,   177.501, 0.482578, 0.628299),
  "SER" => (502, 48.8819, 87.0889, 77.8733,  141.471, 0.561288, 0.550453),
  "HIS" => (155, 42.671,  85.0285, 149.504,  216.355, 0.501843, 0.69101),
  "PRO" => (281, 42.7962, 89.859,   98.9088, 176.266, 0.47626, 0.561134),
  "PHE" => (239, 41.5454, 85.956,  171.084,  238.035, 0.483333, 0.718734),
  "ASP" => (399, 44.0283, 87.96,   118.187, 182.095, 0.500549, 0.649038),
  "ILE" => (338, 40.4473, 86.4206, 136.717, 202.6, 0.468029, 0.674812),
  "TYR" => (299, 42.7139, 85.8969, 187.116, 255.919, 0.497269, 0.731156),
  "ARG" => (230, 45.136,  86.3419, 193.846, 257.924, 0.522759, 0.751562),
  "LEU" => (485, 42.5356, 86.8573, 138.733, 200.877, 0.489719, 0.690638),
  "ALA" => (510, 53.0357, 87.6011, 59.1082, 120.762, 0.605423, 0.489458),
  "MET" => (116, 44.8261, 86.0533, 148.607, 212.591, 0.520911, 0.699026),
  "CYS" => (156, 46.4347, 86.3576, 96.9218, 160.629, 0.537703, 0.603391),
  "GLU" => (341, 45.2395, 87.0802, 142.754, 206.045, 0.519516, 0.692831),
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
                    f_sc = bb_sc_exp()[aa][7]
                    f_bb = bb_sc_exp()[aa][6]
                    GTFEapp = tfe_sc_bb_AutonBolen[aa][i] # apparent free energy
                    f_bb_gly = bb_sc_exp()["GLY"][6]  # = 1
                    tfe_bb = tfe_sc_bb_AutonBolen["BB"][i]
                    γG_val = γG()[cs[i]]              # e.g., +14.47 for urea
                    tfe_sc = (1/f_sc) * (GTFEapp + γG_val + tfe_bb * (f_bb_gly - f_bb))
                    #inv(sc_exp()[aa]) * (
                    #    tfe_sc_bb_AutonBolen[aa][i] + (1 - bb_exp()[aa]) * tfe_gly()[cs[i]] + bb_exp()[aa] * γG()[cs[i]] 
                    #)
                    #(tfe_sc_bb_AutonBolen[aa][i] + γG()[cs[i]] + (1-bb_exp()[aa])*tfe_sc_bb_AutonBolen["BB"][i])/sc_exp()[aa]
                    #tfe_sc_bb_AutonBolen[aa][i] * sc_exp()[aa] 
                    #tfe_sc_bb_AutonBolen[aa][i] + (1-bb_exp0()[aa])*tfe_sc_bb_AutonBolen["BB"][i] + bb_exp0()[aa]*γG()[cs[i]]
                end,
            9)
            for aa in keys(bb_sc_exp())  
        ]...
    )
    data["BB"] = ( 90, 52, 67, 48, 35, 62, -39, 14, 62)
    return data
end

tfe_sc_bb(::Type{MoeserHorinekApp2}) = tfe_sc_bb_MoeserHorinekApp2()

bb_exp0() = Dict{String,Float32}(
    "PHE" => 0.435868,
    "GLN" => 0.429058,
    "ASP" => 0.459705,
    "LYS" => 0.439274,
    "ILE" => 0.350738,
    "TYR" => 0.439274,
    "GLY" => 1.0,
    "ASN" => 0.4563,
    "ARG" => 0.443814,
    "LEU" => 0.400681,
    "TRP" => 0.424518,
    "ALA" => 0.524404,
    "THR" => 0.430193,
    "VAL" => 0.409762,
    "MET" => 0.438138,
    "SER" => 0.499432,
    "GLU" => 0.429058,
    "PRO" => 0.404086,
    "HIS" => 0.45857,
    "CYS" => 0.483541,
)
