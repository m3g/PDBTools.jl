#
# Function to return if an atom is a protein atom from the 
# residue name
#
function isprotein(atom :: AtomType; newres = nothing)
  if atom.resname == newres
    true
  end
  iaa = which_natural_aminoacid(atom)
  if iaa == nothing
    false
  else
    true
  end
end
