function isacidic(atom :: AtomType)
  iaa = which_natural_aminoacid(atom)
  if iaa == nothing
    false
  else
    if natural_aminoacids[iaa].type == "Acidic"
      true
    else
      false
    end
  end
end
