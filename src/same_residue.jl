#
# Function that checks if two atoms belong to the same residue
# without, of course, checking the residue counter
#
function same_residue(atom1 :: Atom, atom2 :: Atom)
  atom1.resnum != atom2.resnum && return false
  atom1.model != atom2.model && return false
  atom1.chain != atom2.chain && return false
  atom1.resname != atom2.resname && return false
  atom1.segname != atom2.segname && return false
  return true
end
