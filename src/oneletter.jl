#
# Function to return a oneletter residue name from the three letter name
#

function oneletter(resname)

  if resname == "ALA" ; return "A" ; end
  if resname == "ARG" ; return "R" ; end
  if resname == "ASN" ; return "N" ; end
  if resname == "ASP" ; return "D" ; end
  if resname == "ASX" ; return "B" ; end
  if resname == "CYS" ; return "C" ; end
  if resname == "GLU" ; return "E" ; end
  if resname == "GLN" ; return "Q" ; end
  if resname == "GLX" ; return "Z" ; end
  if resname == "GLY" ; return "G" ; end
  if resname == "HIS" ; return "H" ; end
  if resname == "HSD" ; return "H" ; end
  if resname == "HSE" ; return "H" ; end
  if resname == "ILE" ; return "I" ; end
  if resname == "LEU" ; return "L" ; end
  if resname == "LYS" ; return "K" ; end
  if resname == "MET" ; return "M" ; end
  if resname == "PHE" ; return "F" ; end
  if resname == "PRO" ; return "P" ; end
  if resname == "SER" ; return "S" ; end
  if resname == "THR" ; return "T" ; end
  if resname == "TRP" ; return "W" ; end
  if resname == "TYR" ; return "Y" ; end
  if resname == "VAL" ; return "V" ; end
  if resname == "HOH" ; return "W" ; end
  if resname == "TIP3P" ; return "W" ; end
  if resname == "TIP4P" ; return "W" ; end
  if resname == "TIP3P" ; return "W" ; end
  if resname == "TIP3" ; return "W" ; end
  if resname == "SPC" ; return "W" ; end

  return "X"
  
end




