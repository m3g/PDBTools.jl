#
# Function that writes a PDB file given the vector of atoms
#

function writePDB( atoms :: Union{Vector{Atom},Vector{MutableAtom}}, filename )

  file = open(filename,"w")
  for atom in atoms
    println(file,write_atom(atom))
  end
  close(file)

end
