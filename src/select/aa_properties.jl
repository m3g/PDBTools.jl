
struct AA
  name::String
  three_letter_code::String
  one_letter_code::String
  type::String
  polar::Bool
  hydrophobic::Bool
  charge::Int
  mono_isotopic_mass::Float64
  mass::Float64
end

const ALA = AA("Alanine"       , "ALA", "A",   "Aliphatic" , false, false,  0,   71.037114,  71.0779 )
const ARG = AA("Arginine"      , "ARG", "R",   "Basic"     , true , false,  1,  156.101111, 156.1857 )
const ASN = AA("Asparagine"    , "ASN", "N",   "Amide"     , true , false,  0,  114.042927, 114.1026 )
const ASP = AA("Aspartic acid" , "ASP", "D",   "Acidic"    , true , false, -1,  115.026943, 115.0874 )
const CYS = AA("Cysteine"      , "CYS", "C",   "Sulfuric"  , false, false,  0,  103.009185,	103.1429 )
const GLN = AA("Glutamine"     , "GLN", "Q",   "Amide"     , true , false,  0,  128.058578,	128.1292 )
const GLU = AA("Glutamic acid" , "GLU", "E",   "Acidic"    , true , false, -1,  129.042593,	129.1140 )
const GLY = AA("Glycine"       , "GLY", "G",   "Aliphatic" , false, false,  0,   57.021464,	 57.0513 )
const HIS = AA("Histidine"     , "HIS", "H",   "Aromatic"  , true , false,  0,  137.058912,	137.1393 )
const HSD = AA("Histidine"     , "HSD", "H",   "Aromatic"  , true , false,  1,  138.067000, 138.1470 )
const HSE = AA("Histidine"     , "HSE", "H",   "Aromatic"  , true , false,  1,  138.067000, 138.1470 )
const ILE = AA("Isoleucine"    , "ILE", "I",   "Aliphatic" , false, true,   0,  113.084064,	113.1576 )
const LEU = AA("Leucine"       , "LEU", "L",   "Aliphatic" , false, true,   0,  113.084064,	113.1576 )
const LYS = AA("Lysine"        , "LYS", "K",   "Basic"     , true , false,  1,  128.094963,	128.1723 )
const MET = AA("Methionine"    , "MET", "M",   "Sulfuric"  , false, false,  0,  131.040485,	131.1961 )
const PHE = AA("Phenylalanine" , "PHE", "F",   "Aromatic"  , false, true,   0,  147.068414,	147.1739 )
const PRO = AA("Proline"       , "PRO", "P",   "Cyclic"    , false, false,  0,   97.052764,	97.1152  )
const SER = AA("Serine"        , "SER", "S",   "Hydroxylic", true , false,  0,   87.032028,	87.07730 )
const THR = AA("Threonine"     , "THR", "T",   "Hydroxylic", true , false,  0,  101.047679,	101.1039 )
const TRP = AA("Tryptophan"    , "TRP", "W",   "Aromatic"  , false, true,   0,  186.079313,	186.2099 )
const TYR = AA("Tyrosine"      , "TYR", "Y",   "Aromatic"  , true , false,  0,  163.063320, 163.1733 )
const VAL = AA("Valine"        , "VAL", "V",   "Aliphatic" , false, true,   0,   99.068414,	 99.1311 )
const natural_aminoacids = [ ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, HSD, 
                             HSE, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, 
                             TYR, VAL ]

"""

```
Sequence
```

Wrapper for strings, or vectors of chars, strings, or residue names, to dispatch on 
functions that operate on amino acid sequences.

# Example

```julia-repl
julia> seq = ["Alanine", "Glutamic acid", "Glycine"];

julia> mass(Sequence(seq))
257.2432

julia> seq = "AEG";

julia> mass(Sequence(seq))
257.2432
```

"""
struct Sequence{T}
  s::T
end

"""

```
mass(s::Sequence)
```

Returns the mass of a sequence of amino acids, given a `Sequence` struct type.

# Examples

```julia-repl
julia> seq = ["Alanine", "Glutamic acid", "Glycine"];

julia> mass(Sequence(seq))
257.2432

julia> seq = "AEG";

julia> mass(Sequence(seq))
257.2432

julia> seq = ["ALA", "GLU", "GLY"];

julia> mass(Sequence(seq))
257.2432
```

"""
function mass(s::Sequence)
  m = 0.
  for aa in s.s
    if length(aa) == 1 
      iaa = findfirst(x -> x.one_letter_code == string(aa), natural_aminoacids)
    elseif length(aa) == 3
      iaa = findfirst(x -> x.three_letter_code == string(aa), natural_aminoacids)
    else
      iaa = findfirst(x -> x.name == string(aa), natural_aminoacids)
    end
    if isnothing(iaa)
        error("Could not residue $aa in natural amino acid list.")
    end
    m += natural_aminoacids[iaa].mass
  end
  return m
end

