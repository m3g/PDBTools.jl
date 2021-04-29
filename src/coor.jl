"""

```
coor(atoms; selection, options) 
```

Returns the coordinates of the atoms. The input may be one atom (type `Atom`), a vector of atoms, or a `Residue`. For vectors or residues, selection strings or functions may be applied. The result may be row-major or column-major. Row-major means that the 3D coordinates will be the three columns of the matrix, thus likely the highest dimension will be the first, if the number of atoms of the selection is greater than 3. 

### Examples

```julia-repl
julia> protein = wget("1LBD");

julia> coor(protein[1])
(45.228, 84.358, 70.638)

julia> coor(protein,"index <= 2")
2×3 Matrix{Float64}:
 45.228  84.358  70.638
 46.08   83.165  70.327

julia> coor(protein,"residue = 1",row_major=false)
3×6 Matrix{Float64}:
 45.228  46.08   45.257  45.823  47.147  46.541
 84.358  83.165  81.872  80.796  82.98   82.639
 70.638  70.327  70.236  69.974  71.413  72.662

julia> coor(protein, only = at -> at.resname == "ALA", row_major=false)
3×110 Matrix{Float64}:
 43.94   43.02   41.996  41.228  42.293  28.854  27.821  26.415  …  17.282  -16.085  -16.377  -17.866  -18.496  -15.888
 81.982  80.825  80.878  79.937  80.676  82.285  81.232  81.715     81.631   84.599   84.019   84.088   83.942   82.583
 70.474  70.455  69.34   69.157  71.812  60.096  60.093  60.408     80.842   53.001   51.738   51.741   52.777   51.706

julia> residues = collect(eachresidue(protein));

julia> coor(residues[1])
6×3 Matrix{Float64}:
 45.228  84.358  70.638
 46.08   83.165  70.327
 45.257  81.872  70.236
 45.823  80.796  69.974
 47.147  82.98   71.413
 46.541  82.639  72.662

```

"""
coor(atom::Atom) = atom.x, atom.y, atom.z

function coor(atoms::AbstractVector{Atom}, selection::String; row_major::Bool=true)
  query = parse_query(selection)
  return coor(atoms,only = atom -> apply_query(query,atom), row_major=row_major)
end

function coor(atoms::AbstractVector{Atom}; only=all, row_major::Bool=true)
  n = 0
  for atom in atoms
    if only(atom)
      n = n + 1
    end
  end 
  if row_major
    x = Matrix{Float64}(undef,n,3)
  else
    x = Matrix{Float64}(undef,3,n)
  end
  i = 0
  for atom in atoms
    if only(atom)
      i = i + 1
      if row_major
        x[i,1:3] .= coor(atom)
      else
        x[1:3,i] .= coor(atom)
      end
    end
  end 
  return x
end

#
# Coordinates of the atoms of a residue/molecule
#
coor(residue::Residue; only=all, row_major::Bool=true) = 
  coor(residue.atoms[residue.range], only=only, row_major=row_major)
coor(residue::Residue, selection::String; row_major::Bool=true) = 
  coor(residue.atoms[residue.range], selection, row_major=row_major)
