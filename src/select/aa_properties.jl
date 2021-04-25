
struct AA
  name::String
  three_letter_code::String
  one_letter_code::String
  type::String
  polar::Bool
  hydrophobic::Bool
  charge::Int
end

ALA = AA("Alanine"       , "ALA", "A",   "Aliphatic" , false, false,  0)
ARG = AA("Arginine"      , "ARG", "R",   "Basic"     , true , false,  1)
ASN = AA("Asparagine"    , "ASN", "N",   "Amide"     , true , false,  0)
ASP = AA("Aspartic acid" , "ASP", "D",   "Acidic"    , true , false, -1)
CYS = AA("Cysteine"      , "CYS", "C",   "Sulfuric"  , false, false,  0)
GLN = AA("Glutamine"     , "GLN", "Q",   "Amide"     , true , false,  0)
GLU = AA("Glutamic acid" , "GLU", "E",   "Acidic"    , true , false, -1)
GLY = AA("Glycine"       , "GLY", "G",   "Aliphatic" , false, false,  0)
HIS = AA("Histidine"     , "HIS", "H",   "Aromatic"  , true , false,  0)
HSD = AA("Histidine"     , "HSD", "H",   "Aromatic"  , true , false,  1)
HSE = AA("Histidine"     , "HSE", "H",   "Aromatic"  , true , false,  1)
ILE = AA("Isoleucine"    , "ILE", "I",   "Aliphatic" , false, true,   0)
LEU = AA("Leucine"       , "LEU", "L",   "Aliphatic" , false, true,   0)
LYS = AA("Lysine"        , "LYS", "K",   "Basic"     , true , false,  1)
MET = AA("Methionine"    , "MET", "M",   "Sulfuric"  , false, false,  0)
PHE = AA("Phenylalanine" , "PHE", "F",   "Aromatic"  , false, true,   0)
PRO = AA("Proline"       , "PRO", "P",   "Cyclic"    , false, false,  0)
SER = AA("Serine"        , "SER", "S",   "Hydroxylic", true , false,  0)
THR = AA("Threonine"     , "THR", "T",   "Hydroxylic", true , false,  0)
TRP = AA("Tryptophan"    , "TRP", "W",   "Aromatic"  , false, true,   0)
TYR = AA("Tyrosine"      , "TYR", "Y",   "Aromatic"  , true , false,  0)
VAL = AA("Valine"        , "VAL", "V",   "Aliphatic" , false, true,   0)

natural_aminoacids = [ ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, HSD, 
                       HSE, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, 
                       TYR, VAL ]



