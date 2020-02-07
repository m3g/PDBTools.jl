#
# Function to return if an atom is a protein atom from the 
# residue name
#

function isprotein(atom :: Atom)

  if atom.resname == "ALA" ; return true ; end
  if atom.resname == "ARG" ; return true ; end
  if atom.resname == "ASN" ; return true ; end
  if atom.resname == "ASP" ; return true ; end
  if atom.resname == "ASX" ; return true ; end
  if atom.resname == "CYS" ; return true ; end
  if atom.resname == "GLU" ; return true ; end
  if atom.resname == "GLN" ; return true ; end
  if atom.resname == "GLX" ; return true ; end
  if atom.resname == "GLY" ; return true ; end
  if atom.resname == "HIS" ; return true ; end
  if atom.resname == "HSD" ; return true ; end
  if atom.resname == "HSE" ; return true ; end
  if atom.resname == "ILE" ; return true ; end
  if atom.resname == "LEU" ; return true ; end
  if atom.resname == "LYS" ; return true ; end
  if atom.resname == "MET" ; return true ; end
  if atom.resname == "PHE" ; return true ; end
  if atom.resname == "PRO" ; return true ; end
  if atom.resname == "SER" ; return true ; end
  if atom.resname == "THR" ; return true ; end
  if atom.resname == "TRP" ; return true ; end
  if atom.resname == "TYR" ; return true ; end
  if atom.resname == "VAL" ; return true ; end

  return false
  
end




