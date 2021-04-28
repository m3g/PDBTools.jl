
struct AA
  name::String
  three_letter_code::String
  one_letter_code::String
  type::String
  polar::Bool
  hydrophobic::Bool
  charge::Int
end

const ALA = AA("Alanine"       , "ALA", "A",   "Aliphatic" , false, false,  0)
const ARG = AA("Arginine"      , "ARG", "R",   "Basic"     , true , false,  1)
const ASN = AA("Asparagine"    , "ASN", "N",   "Amide"     , true , false,  0)
const ASP = AA("Aspartic acid" , "ASP", "D",   "Acidic"    , true , false, -1)
const CYS = AA("Cysteine"      , "CYS", "C",   "Sulfuric"  , false, false,  0)
const GLN = AA("Glutamine"     , "GLN", "Q",   "Amide"     , true , false,  0)
const GLU = AA("Glutamic acid" , "GLU", "E",   "Acidic"    , true , false, -1)
const GLY = AA("Glycine"       , "GLY", "G",   "Aliphatic" , false, false,  0)
const HIS = AA("Histidine"     , "HIS", "H",   "Aromatic"  , true , false,  0)
const HSD = AA("Histidine"     , "HSD", "H",   "Aromatic"  , true , false,  1)
const HSE = AA("Histidine"     , "HSE", "H",   "Aromatic"  , true , false,  1)
const ILE = AA("Isoleucine"    , "ILE", "I",   "Aliphatic" , false, true,   0)
const LEU = AA("Leucine"       , "LEU", "L",   "Aliphatic" , false, true,   0)
const LYS = AA("Lysine"        , "LYS", "K",   "Basic"     , true , false,  1)
const MET = AA("Methionine"    , "MET", "M",   "Sulfuric"  , false, false,  0)
const PHE = AA("Phenylalanine" , "PHE", "F",   "Aromatic"  , false, true,   0)
const PRO = AA("Proline"       , "PRO", "P",   "Cyclic"    , false, false,  0)
const SER = AA("Serine"        , "SER", "S",   "Hydroxylic", true , false,  0)
const THR = AA("Threonine"     , "THR", "T",   "Hydroxylic", true , false,  0)
const TRP = AA("Tryptophan"    , "TRP", "W",   "Aromatic"  , false, true,   0)
const TYR = AA("Tyrosine"      , "TYR", "Y",   "Aromatic"  , true , false,  0)
const VAL = AA("Valine"        , "VAL", "V",   "Aliphatic" , false, true,   0)
const natural_aminoacids = [ ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, HSD, 
                             HSE, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, 
                             TYR, VAL ]



