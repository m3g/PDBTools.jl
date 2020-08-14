#
# Function to return if an atom is an atom from a water molecule
#

function iswater(atom :: AtomType)

  if atom.resname == "HOH"   ; return true ; end
  if atom.resname == "OH2"   ; return true ; end
  if atom.resname == "TIP3"  ; return true ; end
  if atom.resname == "TIP3P" ; return true ; end
  if atom.resname == "TIP4P" ; return true ; end
  if atom.resname == "TIP5P" ; return true ; end
  if atom.resname == "TIP7P" ; return true ; end
  if atom.resname == "SPC"   ; return true ; end
  if atom.resname == "SPCE"  ; return true ; end

  return false
  
end




