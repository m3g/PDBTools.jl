"""

```
getseq(Vector{Atom} or filename; selection, code)
```

Returns the sequence of aminoacids from the vector of atoms or file name. Selections may be applied. Returns both three- and one-letter residue codes. Code defines if the output will be a one-letter, three-letter or full-residue name array.

### Example

```julia-repl
julia> protein = wget("1LBD");

julia> getseq(protein,"residue < 3")
2-element Vector{String}:
 "S"
 "A"

julia> getseq(protein,"residue < 3",code=3)
2-element Vector{String}:
 "Serine"
 "Alanine"

```

"""
function getseq(atoms::AbstractVector{Atom}, selection::String; code::Int=1)
  query = parse_query(selection)
  return getseq(atoms, only = atom -> apply_query(query,atom),code=code)
end

function getseq(atoms::AbstractVector{Atom}; only=all, code::Int=1)
  seq = String[]
  for residue in eachresidue(atoms) 
    # If any atom of this residue is in the selection, add it
    consider = false
    for at in residue
      if isprotein(at) && only(at)
        consider = true
        break
      end
    end
    if consider
      if code == 1
        push!(seq,oneletter(residue.name))
      elseif code == 2
        push!(seq,threeletter(residue.name))
      elseif code == 3
        push!(seq,residuename(residue.name))
      end
    end
  end
  seq
end

# From the file name

function getseq(file::String, selection::String; code::Int=1) 
  atoms = readPDB(file)
  return getseq(atoms, selection, code=code)
end

getseq(file::String; only=all, code::Int=1) = 
   getseq(readPDB(file), only = only, code=code)

