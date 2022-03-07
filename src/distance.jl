"""

```
closest(x,y)
```

Computes the minimum distance between two sets of atoms and returns the indexes of the atoms 
and their distance. Both vector of atoms or vectors of coordinates can be used as input.

### Examples

```julia-repl
julia> model = wget("1BSX");

julia> protein = select(model,"protein");

julia> ligand = select(model,"resname T3");

julia> closest(ligand,protein)
(43, 3684, 2.7775834820937417)

julia> ligand[43]
    4037   O1      T3     B        2      512  -22.568   81.625    3.159 36.59  1.00     1       -      4041

julia> closest(ligand[43],protein)
(1, 3684, 2.7775834820937417)

julia> x = coor(protein)
3994-element Vector{SVector{3, Float64}}:
 [52.884, 24.022, 35.587]
 [52.916, 24.598, 36.993]
 ⋮
 [-46.887, 86.925, 13.235]
 [-47.164, 83.593, 15.25]

julia> closest(ligand,x)
(43, 3684, 2.7775834820937417)

```

"""
function closest(x::AbstractVector{T1}, y::AbstractVector{T2}) where {T1,T2 <: Union{SVector, Atom}}
    imin = -1
    jmin = -1
    dmin = +Inf
    for (i,xatom) in pairs(x) 
        for (j,yatom) in pairs(y)
            d = distance(xatom, yatom)
            if d < dmin
                imin = i
                jmin = j
                dmin = d
            end
        end
    end
    return imin, jmin, dmin 
end
closest(x::Atom, y::AbstractVector{<:Union{SVector, Atom}}) = closest(SVector(x),y)
closest(x::AbstractVector{<:Union{SVector, Atom}}, y::Atom) = closest(x,SVector(y))
closest(x::AbstractVector{<:Union{SVector, Atom}}, y::SVector{3,<:Real}) = closest(x,SVector{1,typeof(y)}(y))
closest(x::SVector{3,<:Real},y::AbstractVector{<:Union{SVector, Atom}}) = closest(SVector{1,typeof(x)}(x),y)
closest(x::Atom, y::Atom) = closest(SVector(x),SVector(y))
closest(x::SVector{3,<:Real}, y::SVector{3,<:Real}) = closest(SVector{1,typeof(x)}(x),SVector{1,typeof(y)}(y))
closest(x::Atom, y::SVector{3,<:Real}) = closest(SVector(x),SVector{1,typeof(y)}(y))
closest(x::SVector{3,<:Real}, y::Atom) = closest(SVector{1,typeof(x)}(x),SVector(y))
closest(x::Residue,y::Residue) = closest(x.atoms[x.range],y.atoms[y.range])

"""

```
distance(x,y)
```

Computes the minimum distance between two sets of atoms, between an atom and a set of atoms, or simply 
the distance between two atoms. The input may be a vector of `Atom`s, or the 
coordinates that are output of the `coor` function. 

### Examples

```julia-repl
julia> model = wget("1BSX");

julia> protein = select(model,"protein");

julia> ligand = select(model,"resname T3");

julia> distance(protein,ligand)
2.7775834820937417

julia> distance(protein[1],ligand[3])
36.453551075306784

julia> distance(coor(ligand),protein)
2.7775834820937417

```

"""
distance(x::SVector, y::SVector) = norm(x - y)
distance(x::Atom, y::Atom) = norm(coor(x) - coor(y))
distance(x::Atom, y::SVector) = norm(coor(x) - y)
distance(x::SVector, y::Atom) = norm(x - coor(y))
distance(x,y) = closest(x,y)[3]

