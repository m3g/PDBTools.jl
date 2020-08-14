function isaliphatic(atom :: AtomType)
  iaa = which_natural_aminoacid(atom)
  if iaa == nothing
    false
  else
    if natural_aminoacids[iaa].type == "Aliphatic"
      true
    else
      false
    end
  end
end
