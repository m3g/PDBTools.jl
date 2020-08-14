#
# Return the coordinates of the atoms
#

function coor( atoms :: Union{Vector{Atom},Vector{MutableAtom}}, selection :: String  )
  query = parse_query(selection)
  return coor(atoms,only = atom -> apply_query(query,atom))
end

function coor( atoms :: Union{Vector{Atom},Vector{MutableAtom}}; only = atom -> true )

  n = 0
  for atom in atoms
    if apply_query(query,atom)
      n = n + 1
    end
  end 
  x = Matrix{Float64}(undef,n,3)
  i = 0
  for atom in atoms
    if only(atom)
      i = i + 1
      x[i,1] = atom.x
      x[i,2] = atom.y
      x[i,3] = atom.z
    end
  end 
  return x
end
