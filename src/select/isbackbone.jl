function isbackbone(atom::Atom)
  iaa = which_natural_aminoacid(atom)
  if iaa == 0
    false
  else
    if atom.name in [ "N", "CA", "C", "O" ]
      true
    else
      false
    end
  end
end
