function isbackbone(atom :: AtomType)
  iaa = which_natural_aminoacid(atom)
  if iaa == nothing
    false
  else
    if atom.name in [ "N", "CA", "C", "O" ]
      true
    else
      false
    end
  end
end
