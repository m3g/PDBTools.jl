#
# Return the coordinates of the atoms
#

struct MaxMinCoords
  xmin :: Vector{Float64}
  xmax :: Vector{Float64}
  xlength :: Vector{Float64}
end

function maxmin( atoms :: AtomVector, selection :: String )
  query = parse_query(selection)
  return maxmin(atoms, only = atom -> apply_query(query,atom))
end

function maxmin( atoms :: AtomVector; only = atom -> true )
  x = coor(atoms; only = only)
  xmin = [ minimum(x[1,:]), minimum(x[2,:]), minimum(x[3,:]) ]
  xmax = [ maximum(x[1,:]), maximum(x[2,:]), maximum(x[3,:]) ]
  xlength = @. xmax - xmin
  return MaxMinCoords(xmin,xmax,xlength)
end

function Base.show( io :: IO, m :: MaxMinCoords )
  println(" ")
  println(" Minimum atom coordinates: xmin = ", m.xmin)
  println(" Maximum atom coordinates: xmax = ", m.xmax)
  println(" Length in each direction: xlength = ", m.xlength)
end

