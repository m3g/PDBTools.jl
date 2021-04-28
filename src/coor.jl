"""

```
coor(atoms; selection, options) 
```

Returns the coordinates of the atoms. The input may be one atom (type `Atom`), a vector of atoms, or a `Residue`. For vectors or residues, selection strings or functions may be applied. The result may be column-based or row based. 

### Examples

```julia-repl
julia> protein = wget("1LBD");

julia> coor(protein[1])
(45.228, 84.358, 70.638)

julia> coor(protein,"residue = 1")
3×6 Matrix{Float64}:
 45.228  46.08   45.257  45.823  47.147  46.541
 84.358  83.165  81.872  80.796  82.98   82.639
 70.638  70.327  70.236  69.974  71.413  72.662

julia> coor(protein, only = at -> at.resname == "ALA" )
3×110 Matrix{Float64}:
 43.94   43.02   41.996  41.228  42.293  28.854  27.821  26.415  …  17.282  -16.085  -16.377  -17.866  -18.496  -15.888
 81.982  80.825  80.878  79.937  80.676  82.285  81.232  81.715     81.631   84.599   84.019   84.088   83.942   82.583
 70.474  70.455  69.34   69.157  71.812  60.096  60.093  60.408     80.842   53.001   51.738   51.741   52.777   51.706

julia> coor(protein,"index <= 2", column_based=true)
3×2 Matrix{Float64}:
 45.228  46.08
 84.358  83.165
 70.638  70.327

julia> residues = collect(eachresidue(protein));

julia> coor(residues[1])
3×6 Matrix{Float64}:
 45.228  46.08   45.257  45.823  47.147  46.541
 84.358  83.165  81.872  80.796  82.98   82.639
 70.638  70.327  70.236  69.974  71.413  72.662

```

"""
coor(atom::Atom) = atom.x, atom.y, atom.z

function coor(atoms::AbstractVector{Atom}, selection::String; column_based::Bool=true)
  query = parse_query(selection)
  return coor(atoms,only = atom -> apply_query(query,atom), column_based=column_based)
end

function coor(atoms::AbstractVector{Atom}; only=all, column_based::Bool=true)
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
        x[1:3,i] .= coor(atom)
      else
        x[i,1:3] .= coor(atom)
      end
    end
  end 
  return x
end

#
# Coordinates of the atoms of a residue/molecule
#
coor(residue::Residue; only=all, column_based::Bool=true) = 
  coor(residue.atoms[residue.range], only=only, column_based=column_based)
coor(residue::Residue, selection::String; column_based::Bool=true) = 
  coor(residue.atoms[residue.range], selection, column_based=column_based)
