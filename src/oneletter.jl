"""
```
oneletter(residue::String)
```

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
    l = length(residue)
    code = uppercase(residue)
    if l > 1
        if l == 3
            ires = findfirst(r -> r.three_letter_code == code, natural_aminoacids)
        else
            ires = findfirst(r -> uppercase(r.name) == code, natural_aminoacids)
        end
        code = (isnothing(ires) ? nothing : natural_aminoacids[ires].one_letter_code)
    end
    return code
end
