#
# Return the coordinates of the atoms
#
coor(atom::Atom) = atom.x, atom.y, atom.z

function coor(atoms::Vector{Atom}, selection::String; column_based::Bool=true)
  query = parse_query(selection)
  return coor(atoms,only = atom -> apply_query(query,atom), column_based=column_based)
end

function coor(atoms::Vector{Atom}; only=all, column_based::Bool=true)
  n = 0
  for atom in atoms
    if only(atom)
      n = n + 1
    end
  end 
  if column_based
    x = Matrix{Float64}(undef,3,n)
  else
    x = Matrix{Float64}(undef,n,3)
  end
  i = 0
  for atom in atoms
    if only(atom)
      i = i + 1
      if column_based 
        x[1:3,i] .= coor(atom)
      else
        x[i,1:3] .= coor(atom)
      end
    end
  end 
  return x
end

# Coordinates of the atoms of a residue/molecule
coor(residue::Residue; column_based::Bool=true) = 
  coor(residue.atoms, only = atom -> (atom.index in residue.range), column_based=column_based)
