"""

```
distance(x,y)
```

Computes the minimum distance between two sets of atoms, between an atom and a set of atoms, or simply the distance between two atoms. The input may be a vector of `Atom`s, or the matrices that are output of the `coor` function. Using the matrices results in greater speed, particularly if the greatest dimension of the matrices is the first one (`xyz_in_cols=true` option, which is the default option).

### Examples

```julia-repl
julia> model = wget("1BSX");

julia> protein = select(model,"protein");

julia> ligand = select(model,"resname T3");

julia> distance(protein,ligand)
2.7775834820937417

julia> distance(protein[1],ligand[3])
36.453551075306784

julia> using BenchmarkTools

julia> @btime distance(\$protein,\$ligand)
  653.590 μs (0 allocations: 0 bytes)
2.7775834820937417

julia> xprot = coor(protein);

julia> xlig = coor(ligand);

julia> @btime distance(\$xprot,\$xlig)
  312.062 μs (0 allocations: 0 bytes)
2.7775834820937417

```

"""
@inline distance(x₁,y₁,z₁,x₂,y₂,z₂) = sqrt((x₂-x₁)^2+(y₂-y₁)^2+(z₂-z₁)^2) 
@inline distance(x,y,z,at::Atom) = distance(x,y,z,at.x,at.y,at.z)
@inline distance(at::Atom,x,y,z) = distance(x,y,z,at)
@inline distance(x::Atom,y::Atom) = distance(x.x,x.y,x.z,y.x,y.y,y.z)

function distance(x::Atom, y::AbstractVector{Atom})
  d = +Inf
  for at in y
    d = min(d,distance(x,at))
  end
  d
end

function distance(x::AbstractVector{Atom}, y::AbstractVector{Atom})
  d = +Inf
  for atx in x
    d = min(d,distance(atx,y))
  end
  d
end

distance(r1::Residue,r2::Residue) = 
  distance(@view(r1.atoms[r1.range]),@view(r2.atoms[r2.range]))

using LoopVectorization

distance(at::Atom,A::Matrix{T}; xyz_in_cols=true) where T =
  distance(at.x,at.y,at.z,A,xyz_in_cols=xyz_in_cols)

function distance(x,y,z,A::Matrix{T}; xyz_in_cols=true) where T
  n, m = size(A)
  d = +Inf
  if xyz_in_cols
    @assert m == 3 "Number of columns of coordinates matrix must be 3 when xyz_in_cols=true"
    @avx for i in 1:n
      d = min(d,distance(x,y,z,A[i,1],A[i,2],A[i,3]))
    end
  else
    @assert n == 3 "Number of rows of coordinates matrix must be 3 when xyz_in_cols=false"
    @avx for j in 1:m
      d = min(d,distance(x,y,z,A[1,j],A[2,j],A[3,j]))
    end
  end
  d
end

function distance(A::Matrix{T},B::Matrix{S}; xyz_in_cols=true) where {T,S}
  nA, mA = size(A)
  nB, mB = size(B)
  d = +Inf
  if xyz_in_cols
    @assert (mA == 3 && mB == 3) "Matrices must have 3 columns if xyz_in_cols=true"
    if nA > nB
      for i in 1:nB
        d = min(d,distance(B[i,1],B[i,2],B[i,3],A,xyz_in_cols=true))
      end
    else
      for i in 1:nA
        d = min(d,distance(A[i,1],A[i,2],A[i,3],B,xyz_in_cols=true))
      end
    end
  else
    @assert (nA == 3 && nB == 3) "Matrices must have 3 rows if xyz_in_cols=false" 
    if mA > mB
      for j in 1:mB
        d = min(d,distance(B[1,j],B[2,j],B[3,j],A,xyz_in_cols=false))
      end
    else
      for j in 1:mA
        d = min(d,distance(A[1,j],A[2,j],A[3,j],B,xyz_in_cols=false))
      end
    end
  end
  d
end

