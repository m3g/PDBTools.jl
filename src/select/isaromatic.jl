function isaromatic(atom :: AtomType)
  iaa = which_natural_aminoacid(atom)
  if iaa == nothing
    false
  else
    if natural_aminoacids[iaa].type == "Aromatic"
      true
    else
      false
    end
  end
end
