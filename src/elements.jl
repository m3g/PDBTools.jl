#
# List of elements with properties
#
struct Element
 atomic_number::Int
 element::String
 pdb_name::String
 name::String
 mass::Float64
end
const elements = [
             Element( 0      ,  "X"  ,  "X"   ,   "NotFound"    ,  0.00000   )
             Element( 1      ,  "H"  ,  "H"   ,   "Hydrogen"    ,  1.00797   )
             Element( 2      ,  "He" ,  "HEL" ,   "Helium"      ,  4.00260   )
             Element( 3      ,  "Li" ,  "LI"  ,   "Lithium"     ,  6.941     )
             Element( 4      ,  "Be" ,  "BE"  ,   "Beryllium"   ,  9.01218   )
             Element( 5      ,  "B"  ,  "B"   ,   "Boron"       ,  10.81     )
             Element( 6      ,  "C"  ,  "C"   ,   "Carbon"      ,  12.011    )
             Element( 7      ,  "N"  ,  "N"   ,   "Nitrogen"    ,  14.0067   )
             Element( 8      ,  "O"  ,  "O"   ,   "Oxygen"      ,  15.9994   )
             Element( 9      ,  "F"  ,  "F"   ,   "Fluorine"    ,  18.998403 )
             Element( 10     ,  "Ne" ,  "NEO" ,   "Neon"        ,  20.179    )
             Element( 11     ,  "Na" ,  "NA"  ,   "Sodium"      ,  22.98977  )
             Element( 11     ,  "Na" ,  "SOD" ,   "Sodium"      ,  22.98977  )
             Element( 12     ,  "Mg" ,  "MG"  ,   "Magnesium"   ,  24.305    )
             Element( 13     ,  "Al" ,  "AL"  ,   "Aluminum"    ,  26.98154  )
             Element( 14     ,  "Si" ,  "SI"  ,   "Silicon"     ,  28.0855   )
             Element( 15     ,  "P"  ,  "P"   ,   "Phosphorus"  ,  30.97376  )
             Element( 16     ,  "S"  ,  "S"   ,   "Sulfur"      ,  32.06     )
             Element( 17     ,  "Cl" ,  "CL"  ,   "Chlorine"    ,  35.453    )
             Element( 17     ,  "Cl" ,  "CLA" ,   "Chlorine"    ,  35.453    )
             Element( 18     ,  "Ar" ,  "AR"  ,   "Argon"       ,  39.948    )
             Element( 19     ,  "K"  ,  "K"   ,   "Potassium"   ,  39.0983   )
             Element( 19     ,  "K"  ,  "POT" ,   "Potassium"   ,  39.0983   )
             Element( 20     ,  "Ca" ,  "CAL" ,   "Calcium"     ,  40.08     )
             Element( 21     ,  "Sc" ,  "SC"  ,   "Scandium"    ,  44.9559   )
             Element( 22     ,  "Ti" ,  "TI"  ,   "Titanium"    ,  47.90     )
             Element( 23     ,  "V"  ,  "V"   ,   "Vanadium"    ,  50.9415   )
             Element( 24     ,  "Cr" ,  "CR"  ,   "Chromium"    ,  51.996    )
             Element( 25     ,  "Mn" ,  "MN"  ,   "Manganese"   ,  54.9380   )
             Element( 26     ,  "Fe" ,  "FE"  ,   "Iron"        ,  55.847    )
             Element( 27     ,  "Co" ,  "CO"  ,   "Cobalt"      ,  58.9332   )
             Element( 28     ,  "Ni" ,  "NI"  ,   "Nickel"      ,  58.70     )
             Element( 29     ,  "Cu" ,  "CU"  ,   "Copper"      ,  63.546    )
             Element( 30     ,  "Zn" ,  "ZN"  ,   "Zinc"        ,  65.38     )
             Element( 31     ,  "Ga" ,  "GA"  ,   "Gallium"     ,  69.72     )
             Element( 32     ,  "Ge" ,  "GE"  ,   "Germanium"   ,  72.59     )
             Element( 33     ,  "As" ,  "AS"  ,   "Arsenic"     ,  74.9216   )
             Element( 34     ,  "Se" ,  "SE"  ,   "Selenium"    ,  78.96     )
             Element( 35     ,  "Br" ,  "BR"  ,   "Bromine"     ,  79.904    )
             Element( 36     ,  "Kr" ,  "KR"  ,   "Krypton"     ,  83.80     )
             Element( 37     ,  "Rb" ,  "RB"  ,   "Rubidium"    ,  85.4678   )
             Element( 38     ,  "Sr" ,  "SR"  ,   "Strontium"   ,  87.62     )
             Element( 39     ,  "Y"  ,  "Y"   ,   "Yttrium"     ,  88.9059   )
             Element( 40     ,  "Zr" ,  "ZR"  ,   "Zirconium"   ,  91.22     )
             Element( 41     ,  "Nb" ,  "NB"  ,   "Niobium"     ,  92.9064   )
             Element( 42     ,  "Mo" ,  "MO"  ,   "Molybdenum"  ,  95.94     )
             Element( 43     ,  "Tc" ,  "TC"  ,   "Technetium"  ,  98.       )
             Element( 44     ,  "Ru" ,  "RU"  ,   "Ruthenium"   ,  101.07    )
             Element( 45     ,  "Rh" ,  "RH"  ,   "Rhodium"     ,  102.9055  )
             Element( 46     ,  "Pd" ,  "PD"  ,   "Palladium"   ,  106.4     )
             Element( 47     ,  "Ag" ,  "AG"  ,   "Silver"      ,  107.868   )
             Element( 48     ,  "Cd" ,  "CAD" ,   "Cadmium"     ,  112.41    )
             Element( 49     ,  "In" ,  "IN"  ,   "Indium"      ,  114.82    )
             Element( 50     ,  "Sn" ,  "SN"  ,   "Tin"         ,  118.69    )
             Element( 51     ,  "Sb" ,  "SB"  ,   "Antimony"    ,  121.75    )
             Element( 52     ,  "Te" ,  "TE"  ,   "Tellurium"   ,  127.60    )
             Element( 53     ,  "I"  ,  "I"   ,   "Iodine"      ,  126.9045  )
             Element( 54     ,  "Xe" ,  "XE"  ,   "Xenon"       ,  131.30    )
             Element( 55     ,  "Cs" ,  "CES" ,   "Cesium"      ,  132.9054  )
             Element( 56     ,  "Ba" ,  "BA"  ,   "Barium"      ,  137.33    )
             Element( 57     ,  "La" ,  "LA"  ,   "Lanthanum"   ,  138.9055  )
             Element( 58     ,  "Ce" ,  "CER" ,   "Cerium"      ,  140.12    )
             Element( 59     ,  "Pr" ,  "PR"  ,   "Praseodymium",  140.9077  )
             Element( 60     ,  "Nd" ,  "ND"  ,   "Neodymium"   ,  144.24    )
             Element( 61     ,  "Pm" ,  "PM"  ,   "Promethium"  ,  145       )
             Element( 62     ,  "Sm" ,  "SM"  ,   "Samarium"    ,  150.4     )
             Element( 63     ,  "Eu" ,  "EU"  ,   "Europium"    ,  151.96    )
             Element( 64     ,  "Gd" ,  "GD"  ,   "Gadolinium"  ,  157.25    )
             Element( 65     ,  "Tb" ,  "TB"  ,   "Terbium"     ,  158.9254  )
             Element( 66     ,  "Dy" ,  "DY"  ,   "Dysprosium"  ,  162.50    )
             Element( 67     ,  "Ho" ,  "HLM" ,   "Holmium"     ,  164.9304  )
             Element( 68     ,  "Er" ,  "ER"  ,   "Erbium"      ,  167.26    )
             Element( 69     ,  "Tm" ,  "TM"  ,   "Thulium"     ,  168.9342  )
             Element( 70     ,  "Yb" ,  "YB"  ,   "Ytterbium"   ,  173.04    )
             Element( 71     ,  "Lu" ,  "LU"  ,   "Lutetium"    ,  174.967   )
             Element( 72     ,  "Hf" ,  "HAF" ,   "Hafnium"     ,  178.49    )
             Element( 73     ,  "Ta" ,  "TA"  ,   "Tantalum"    ,  180.9479  )
             Element( 74     ,  "W"  ,  "W"   ,   "Tungsten"    ,  183.85    )
             Element( 75     ,  "Re" ,  "RE"  ,   "Rhenium"     ,  186.207   )
             Element( 76     ,  "Os" ,  "OS"  ,   "Osmium"      ,  190.2     )
             Element( 77     ,  "Ir" ,  "IR"  ,   "Iridium"     ,  192.22    )
             Element( 78     ,  "Pt" ,  "PT"  ,   "Platinum"    ,  195.09    )
             Element( 79     ,  "Au" ,  "AU"  ,   "Gold"        ,  196.9665  )
             Element( 80     ,  "Hg" ,  "MER" ,   "Mercury"     ,  200.59    )
             Element( 81     ,  "Tl" ,  "TL"  ,   "Thallium"    ,  204.37    )
             Element( 82     ,  "Pb" ,  "PB"  ,   "Lead"        ,  207.2     )
             Element( 83     ,  "Bi" ,  "BI"  ,   "Bismuth"     ,  208.9804  )
             Element( 84     ,  "Po" ,  "PO"  ,   "Polonium"    ,  209       )
             Element( 85     ,  "At" ,  "AT"  ,   "Astatine"    ,  210       )
             Element( 86     ,  "Rn" ,  "RN"  ,   "Radon"       ,  222       )
             Element( 87     ,  "Fr" ,  "FR"  ,   "Francium"    ,  223       )
             Element( 88     ,  "Ra" ,  "RA"  ,   "Radium"      ,  226.0254  )
             Element( 89     ,  "Ac" ,  "AC"  ,   "Actinium"    ,  227.0278  )
             Element( 90     ,  "Th" ,  "TH"  ,   "Thorium"     ,  232.0381  )
             Element( 91     ,  "Pa" ,  "PA"  ,   "Protactinium",  231.0359  )
             Element( 92     ,  "U"  ,  "U"   ,   "Uranium"     ,  238.029   )
           ]                        


#
# Retrive index of element in elements list from name. Returns 1 (element "X" of list) if not found
#
function element_index(name::String)
  # Try to find if there is any exact match
  i = findfirst(el -> (el.element == name || 
                       el.name == name ||
                       el.pdb_name == name ), elements)
  (i != nothing) && return i
  #
  # Try to find by PDB name, note that NT2 and 2NT2 must match N, for example
  # 
  # First, check if the first char is a number
  i0 = try 
    parse(Int,name[1:1])
    2
  catch
    1
  end
  # Now check if the start of this name matches some PDB code
  for (iel,el) in pairs(elements)
    lpdb = length(el.pdb_name)
    # If the name is shorter than the PDB name, no match
    (length(name) - i0 + 1) < lpdb && continue
    if name[i0:i0+lpdb-1] == el.pdb_name
      return iel
    end
  end
  return 1
end



