function ishydrophobic(atom :: Union{Atom,MutableAtom})
  iaa = which_natural_aminoacid(atom)
  if iaa == nothing
    false
  else
    if natural_aminoacids[iaa].hydrophobic == true
      true
    else
      false
    end
  end
end