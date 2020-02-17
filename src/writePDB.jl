#
# Function that writes a PDB file given the vector of atoms, with
# optional definition of a selection to be print
#

function writePDB( atoms :: Union{Vector{Atom},Vector{MutableAtom}}, filename; sel = Nothing )

  file = open(filename,"w")
  if sel == Nothing
    for atom in atoms
      println(file,write_atom(atom))
    end
  else
    for i in 1:length(sel)
      println(file,write_atom(atoms[sel[i]]))
    end
  end
  close(file)

end
