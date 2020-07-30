#
# Function to return if an atom is a protein atom from the 
# residue name
#

function isprotein(atom :: Union{Atom,MutableAtom}; newres = Nothing)

  # To take into account alternate conformations, such as "AGLY", "BGLY", etc.
  l = length(atom.resname)
  name = atom.resname[l-2:l]

  if name == "ALA" ; return true ; end
  if name == "ARG" ; return true ; end
  if name == "ASN" ; return true ; end
  if name == "ASP" ; return true ; end
  if name == "ASX" ; return true ; end
  if name == "CYS" ; return true ; end
  if name == "GLU" ; return true ; end
  if name == "GLN" ; return true ; end
  if name == "GLX" ; return true ; end
  if name == "GLY" ; return true ; end
  if name == "HIS" ; return true ; end
  if name == "HSD" ; return true ; end
  if name == "HSE" ; return true ; end
  if name == "ILE" ; return true ; end
  if name == "LEU" ; return true ; end
  if name == "LYS" ; return true ; end
  if name == "MET" ; return true ; end
  if name == "PHE" ; return true ; end
  if name == "PRO" ; return true ; end
  if name == "SER" ; return true ; end
  if name == "THR" ; return true ; end
  if name == "TRP" ; return true ; end
  if name == "TYR" ; return true ; end
  if name == "VAL" ; return true ; end

  # If this is a new residue name indicated by the user
  if newres != Nothing
    if typeof(newres) == String 
      if name == newres ; return true ; end 
    elseif typeof(newres) <: Array{String} 
       if name in newres ; return true ; end
    end
  end

  return false
  
end




