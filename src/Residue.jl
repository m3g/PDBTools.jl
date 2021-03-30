"""

```
Residue(atoms::AbstractVector{Atom}, range::UnitRange{Int})
```

Residue data structure. It contains two fields: `atoms` which is a vector of
`Atom` elements, and `range`, which indicates which atoms of the `atoms` vector
compose the residue.

The purpose of this structure is the create an iterator over the residues of a
vector of atoms.

### Example

```julia-repl

julia> atoms = readPDB("../test/structure.pdb","protein");

julia> for res in eachresidue(atoms)
         name(res) == "ALA" && println(res[1])
       end
   index name resname chain   resnum  residue        x        y        z     b occup model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1

   index name resname chain   resnum  residue        x        y        z     b occup model segname index_pdb
     232    N     ALA     A       19       19   -6.690   -3.943    3.356  0.00  1.00     1    PROT       232

```

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
struct EachResidue
  atoms::AbstractVector{Atom}
end
eachresidue(atoms::AbstractVector{Atom}) = EachResidue(atoms)

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

