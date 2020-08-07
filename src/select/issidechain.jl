function issidechain(atom :: Union{Atom,MutableAtom})
  iaa = which_natural_aminoacid(atom)
  if iaa == nothing
    false
  else
    if ! (atom.name in [ "N", "CA", "C", "O", "HN", "HA", "HT1", "HT2", "HT3" ])
      true
    else
      false
    end
  end
end
