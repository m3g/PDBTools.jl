#
# Return the coordinates of the atoms
#

function xAll( atoms :: Vector{Atom} )

  n = length(atoms)
  x = Matrix{Float64}(undef,n,3)
  i = 0
  for atom in atoms
    i = i + 1
    x[i,1] = atom.x
    x[i,2] = atom.y
    x[i,3] = atom.z
  end 

  return x

end
