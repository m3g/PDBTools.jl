"""

```
Residue(atoms::AbstractVector{Atom}, range::UnitRange{Int})
```

Residue data structure. It contains two fields: `atoms` which is a vector of
`Atom` elements, and `range`, which indicates which atoms of the `atoms` vector
compose the residue.

The purpose of this structure is the create an iterator over the residues of a
vector of atoms.


"""
struct Residue{T<:AbstractVector{Atom}}
  atoms::T
  range::UnitRange{Int}
end
name(residue::Residue) = residue[1].resname
range(residue::Residue) = residue.range

#
# Structure and function to define the eachresidue iterator
#
struct EachResidue{T<:AbstractVector{Atom}}
  atoms::T
end

"""

`eachresidue(atoms::AbstractVector{Atom})`

Iterator for the residues (or molecules) of a selection. 

### Example

```julia-repl
julia> atoms = wget("1LBD");

julia> length(eachresidue(atoms))
238

julia> for res in eachresidue(atoms)
         println(res)
       end
 Residue of name SER with 6 atoms.
   index name resname chain   resnum  residue        x        y        z     b occup model segname index_pdb
       1    N     SER     A      225        1   45.228   84.358   70.638 67.05  1.00     1       -         1
       2   CA     SER     A      225        1   46.080   83.165   70.327 68.73  1.00     1       -         2
       3    C     SER     A      225        1   45.257   81.872   70.236 67.90  1.00     1       -         3
       4    O     SER     A      225        1   45.823   80.796   69.974 64.85  1.00     1       -         4
       5   CB     SER     A      225        1   47.147   82.980   71.413 70.79  1.00     1       -         5
       6   OG     SER     A      225        1   46.541   82.639   72.662 73.55  1.00     1       -         6

 Residue of name ALA with 5 atoms.
   index name resname chain   resnum  residue        x        y        z     b occup model segname index_pdb
       7    N     ALA     A      226        2   43.940   81.982   70.474 67.09  1.00     1       -         7
       8   CA     ALA     A      226        2   43.020   80.825   70.455 63.69  1.00     1       -         8
       9    C     ALA     A      226        2   41.996   80.878   69.340 59.69  1.00     1       -         9
                                                      ...

```

"""
eachresidue(atoms::AbstractVector{Atom}) = EachResidue(atoms)

# Collect residues default constructor
Base.collect(r::EachResidue) = collect(Residue,r)

#
# Iterate over the resiudes
#
function Base.iterate(residues::EachResidue, state=1)
  r0 = state
  r0 > length(residues.atoms) && return nothing
  residue0 = residues.atoms[r0].residue
  r1 = r0
  while r1 <= length(residues.atoms)
    if residues.atoms[r1].residue != residue0 
      return (Residue(residues.atoms,r0:r1-1),r1)
    end
    r1 += 1
  end
  return (Residue(residues.atoms,r0:r1-1),r1)
end

#
# Iterate over atoms of one residue
#
function Base.iterate(residue::Residue, state=1)
  i1 = residue.range[begin] + state - 1
  if i1 <= residue.range[end]
    return (residue.atoms[i1],state+1)
  else
    return nothing
  end
end

#
# Length of the eachresidue iterator (number of residues)
#
function Base.length(residues::EachResidue)
  n = 0
  for residue in residues
    n += 1
  end
  n
end

import Base.getindex
function getindex(residue::Residue,i) 
  @assert i > 0 "Index must be in 1:$(residue.range[end]-residue.range[begin])"
  @assert (i <= length(residue.range)) "Residue has $(residue.range[end]-residue.range[begin]+1) atoms."
  i = residue.range[begin] + i - 1
  residue.atoms[i]
end

function Base.show(io::IO, residue::Residue)
  natoms = residue.range[end]-residue.range[begin]+1
  println(" Residue of name $(name(residue)) with $natoms atoms.")
  print_short_atom_list(@view residue.atoms[residue.range])
end

function Base.show(io::IO, residues::EachResidue)
  println(" Iterator with $(length(residues)) residues.")
end

function Base.show( io :: IO,::MIME"text/plain", residues::AbstractVector{Residue} )
  println("   Array{Residue,1} with $(length(residues)) residues.")
end



