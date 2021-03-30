
"""

```
residues(atoms::AbstractVector{Atom})
```

Returns an iterable object along the residues of the `atoms` vector.

### Example

```julia-repl
julia> pep = readPDB("peptide.pdb");

julia> for residue in resiudes(peptide)
         println(residue.name)
       end

```

"""
