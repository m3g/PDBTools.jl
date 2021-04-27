"""
```
threeletter(residue::String) 
```

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
  l = length(residue)
  code = uppercase(residue)
  if l != 3
    if l == 1
      ires = findfirst(r->r.one_letter_code == code, natural_aminoacids)
    else
      ires = findfirst(r->uppercase(r.name) == code, natural_aminoacids)
    end
    code = (ires == nothing ? nothing : natural_aminoacids[ires].three_letter_code)  
  end
  return code
end
