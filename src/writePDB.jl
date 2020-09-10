#
# Function that writes a PDB file given the vector of atoms, with
# optional definition of a selection to be print
#

function writePDB( atoms :: Vector{<:AtomType}, filename, selection)
  query = parse_query(selection)
  writePDB(atoms,filename,only=atom->apply_query(query,atom))
end

function writePDB( atoms :: Vector{<:AtomType}, filename; only = atom -> true)
  file = open(filename,"w")
  for atom in atoms
    if only(atom)
      println(file,write_atom(atom))
    end
  end
  close(file)
end

