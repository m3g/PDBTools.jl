function isacidic(atom :: Atom)
  iaa = which_natural_aminoacid(atom)
  if iaa == 0
    false
  else
    if natural_aminoacids[iaa].type == "Acidic"
      true
    else
      false
    end
  end
end
