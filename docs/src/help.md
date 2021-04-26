# Help entries

These entries can be accessed from the Julia REPL by typing `?`, for example,
```julia-repl
julia> ? mass
search: mass mapslices MathConstants makedocs set_zero_subnormals get_zero_subnormals mutable struct

  mass(name::String or atom::Atom or Vector{Atom})

  Returns the mass of an atom given its name, or Atom structure, or the total mass of a vector of Atoms.

  Example
  –––––––––

  julia> atoms = [ Atom(name="NT3"), Atom(name="CA") ];
  
  julia> mass(atoms[1])
  14.0067
  
  julia> mass("CA")
  12.011
  
  julia> mass(atoms)
  26.017699999999998

```

```@autodocs
Modules=[PDBTools]
```
