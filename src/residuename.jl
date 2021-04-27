"""
```
residuename(residue::String)
```

Function to return the long residue name from other residue codes. The function is case-insensitive.

### Examples

```julia-repl
julia> residuename("A")
"Alanine"

julia> residuename("Glu")
"Glutamic Acid"

```

"""
function residuename(residue::String)
  l = length(residue)
  code = uppercase(residue)
  if l == 1
    ires = findfirst(r->r.one_letter_code == code, natural_aminoacids)
  elseif l == 3
    ires = findfirst(r->r.three_letter_code == code, natural_aminoacids)
  else
    ires = findfirst(r->uppercase(r.name) == code, natural_aminoacids)
  end
  code = (ires == nothing ? nothing : natural_aminoacids[ires].name)
  return code
end

