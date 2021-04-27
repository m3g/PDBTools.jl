"""
```
Formula::DataType
```

Formula data type. Contains the number of atoms of each type in a vector of tuples.

### Example

```julia-repl
julia> atoms = wget("1LBD","protein and residue 1");

julia> f = formula(atoms)
C₃N₁O₂


julia> f[1]
("C", 3)

julia> f[2]
("N", 1)

julia> f[3]
("O", 2)

```

"""
struct Formula
  formula::Vector{Tuple{String,Int}}
end
Base.getindex(f::Formula,i) = f.formula[i]
"""

```
formula(atoms::AbstractVector{Atom})
```

Returns the molecular formula of the current selection. 

### Example

```julia-repl
julia> first_residue = wget("1LBD","protein and residue 1")
   Array{Atoms,1} with 6 atoms with fields:
   index name resname chain   resnum  residue        x        y        z     b occup model segname index_pdb
       1    N     SER     A      225        1   45.228   84.358   70.638 67.05  1.00     1       -         1
       2   CA     SER     A      225        1   46.080   83.165   70.327 68.73  1.00     1       -         2
       3    C     SER     A      225        1   45.257   81.872   70.236 67.90  1.00     1       -         3
       4    O     SER     A      225        1   45.823   80.796   69.974 64.85  1.00     1       -         4
       5   CB     SER     A      225        1   47.147   82.980   71.413 70.79  1.00     1       -         5
       6   OG     SER     A      225        1   46.541   82.639   72.662 73.55  1.00     1       -         6


julia> formula(first_residue)
C₃N₁O₂

```

"""
function formula(atoms::AbstractVector{Atom})
  f = Formula(Tuple{String,Int}[])
  for at in atoms
    el = element(at)
    i = findfirst(elm -> elm[1] == el, f.formula)
    if i == nothing
      push!(f.formula,(el,1))
    else
      f.formula[i] = (el,f.formula[i][2]+1)
    end
  end
  # Sort by atomic number
  sort!(f.formula, by = x -> element_index(x[1]))
  return f
end
function Base.show(io::IO, f::Formula)
  s = ""
  for el in f.formula
    s *= el[1]*format(el[2])
  end
  for sub in sub_int
    s = replace(s,sub)
  end
  println(s)
end

"""
```
stoichiometry(atoms::AbstractVector{Atom})
```

Returns the stoichiometry of atom selection in a `Formula` structure. 

### Example

```julia-repl
julia> stoichiometry(select(atoms,"water"))
H₂O₁

```

"""
function stoichiometry(atoms::AbstractVector{Atom})
  f = formula(atoms)
  d = gcd( (x[2] for x in f.formula)... ) 
  for (i,p) in pairs(f.formula)
    f.formula[i] = (p[1],p[2]÷d)
  end
  f
end

const sub_int = ( "0" => "₀" ,
                  "1" => "₁" ,
                  "2" => "₂" ,
                  "3" => "₃" ,
                  "4" => "₄" ,
                  "5" => "₅" ,
                  "6" => "₆" ,
                  "7" => "₇" ,
                  "8" => "₈" ,
                  "9" => "₉" )

