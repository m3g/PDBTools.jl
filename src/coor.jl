#
# Return the coordinates of the atoms
#

function coor( atoms :: AtomVector, selection :: String; column_based :: Bool = true  )
  query = parse_query(selection)
  return coor(atoms,only = atom -> apply_query(query,atom), column_based = column_based)
end

function coor( atoms :: AtomVector; only = atom -> true, column_based :: Bool = true )
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
        x[1,i] = atom.x
        x[2,i] = atom.y
        x[3,i] = atom.z
      else
        x[i,1] = atom.x
        x[i,2] = atom.y
        x[i,3] = atom.z
      end
    end
  end 
  return x
end
