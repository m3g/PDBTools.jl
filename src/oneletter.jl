"""

`oneletter(residue::String)` 

Function to return a one-letter residue code from the three letter code or residue name. The function is case-insensitive.

### Examples

```julia-repl
julia> oneletter("ALA")
"A"

julia> oneletter("Glutamic acid")
"E"

```

"""
function oneletter(residue::String)
  if length(residue) == 3
    ires = findfirst(r->lowercase(r.three_letter_code) == lowercase(residue), natural_aminoacids)
  else
    ires = findfirst(r->lowercase(r.name) == lowercase(residue), natural_aminoacids)
  end
  if ires == nothing
    return nothing
  else
    return natural_aminoacids[ires].one_letter_code
  end
end

