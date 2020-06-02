#
# Return the coordinates of the atoms with name "CB", except
# for GLY, in which case we return the coordinates of the "CA" atom
#

function xCB( atoms :: Vector{Atom} )

  ngly = count( atom -> atom.resname == "GLY", atoms )
  n = count( atom -> atom.name == "CB", atoms ) + ngly
  x = Matrix{Float64}(undef,n,3)
  i = 0
  for atom in atoms
    if atom.resname == "GLY" 
      if atom.name == "CA"
        i = i + 1
        x[i,1] = atom.x
        x[i,2] = atom.y
        x[i,3] = atom.z
      end
    else
      if atom.name == "CB"
        i = i + 1
        x[i,1] = atom.x
        x[i,2] = atom.y
        x[i,3] = atom.z
      end
    end
  end 
  return x

end
