"""

`threeletter(residue::String)` 

Function to return a three-letter residue code from the one-letter code or residue name. The function is case-insensitive.

### Examples

```julia-repl
julia> threeletter("A")
"ALA"

julia> threeletter("Aspartic acid")
"ASP"

```

"""
function threeletter(residue::String)
  if length(residue) == 1
    ires = findfirst(r->lowercase(r.one_letter_code) == lowercase(residue), natural_aminoacids)
  else
    ires = findfirst(r->lowercase(r.name) == lowercase(residue), natural_aminoacids)
  end
  if ires == nothing
    return nothing
  else
    return natural_aminoacids[ires].three_letter_code
  end
end

