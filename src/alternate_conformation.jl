# This function will return 'nothing' for alternate (BGLY, CARG, etc) conformations
# of protein residues, which will then be ignored. Only "A" conformations are kept.
# Other non-protein residues are always kept.

function alternate_conformation(atom::Atom)
  if which_natural_aminoacid(atom) != 0
    if length(atom.resname) == 4
      if atom.resname[1:1] == "A"
        return atom.resname[2:4]
      else 
        return nothing
      end
    end
  end
  return atom.resname
end

