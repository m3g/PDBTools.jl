#
# Return the coordinates of the atoms with name "Name"
#

function xName( atoms :: Vector{Atom}, name :: String )

  n = count( atom -> atom.name == name, atoms )
  x = Matrix{Float64}(undef,n,3)
  i = 0
  for atom in atoms
    if atom.name == name
      i = i + 1
      x[i,1] = atom.x
      x[i,2] = atom.y
      x[i,3] = atom.z
    end
  end 

  return x

end
