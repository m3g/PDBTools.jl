#
# Return the coordinates of the atoms
#

struct MaxMinCoords
  xmin::Vector{Float64}
  xmax::Vector{Float64}
  xlength::Vector{Float64}
end

"""

```
maxmin(atoms::Vector{Atom}; selection)
```

Returns the maximum and minimum coordinates of an atom vector, and the length (maximum minus minimum) in each direction. 

### Example

```julia-repl
julia> protein = wget("1LBD");

julia> maxmin(protein)
 
 Minimum atom coordinates: xmin = [-29.301, 57.178, 45.668]
 Maximum atom coordinates: xmax = [47.147, 99.383, 86.886]
 Length in each direction: xlength = [76.448, 42.205, 41.217999999999996]

```



"""
function maxmin(atoms::AbstractVector{Atom}, selection::String)
  query = parse_query(selection)
  return maxmin(atoms, only = atom -> apply_query(query,atom))
end

function maxmin(atoms::AbstractVector{Atom}; only=all)
  x = coor(atoms; only = only)
  xmin = [ minimum(x[1,:]), minimum(x[2,:]), minimum(x[3,:]) ]
  xmax = [ maximum(x[1,:]), maximum(x[2,:]), maximum(x[3,:]) ]
  xlength = @. xmax - xmin
  return MaxMinCoords(xmin,xmax,xlength)
end

function Base.show(io::IO, m::MaxMinCoords)
  println(" ")
  println(" Minimum atom coordinates: xmin = ", m.xmin)
  println(" Maximum atom coordinates: xmax = ", m.xmax)
  println(" Length in each direction: xlength = ", m.xlength)
end

