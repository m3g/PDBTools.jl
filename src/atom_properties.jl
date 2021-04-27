"""
```
atomic_number(name::String or atom::Atom)
```

Returns the atomic number of an atom given its name, or `Atom` structure.

### Example

```julia-repl
julia> at = Atom(name="NT3");

julia> atomic_number(at)
7

julia> atomic_number("CA")
6

```

"""
atomic_number(name::String) = elements[element_index(name)].atomic_number
atomic_number(atom::Atom) = atomic_number(atom.name)

"""
```
element_name(name::String or atom::Atom)
```

Returns the element name of an atom given its name, or `Atom` structure.

### Example

```julia-repl
julia> at = Atom(name="NT3");

julia> element_name(at)
"Nitrogen"

julia> element_name("NT3")
"Nitrogen"

julia> element_name("CA")
"Carbon"

```

"""
element(name::String) = elements[element_index(name)].element
element(atom::Atom) = element(atom.name)
element_name(name::String) = elements[element_index(name)].name  
element_name(atom::Atom) = element_name(atom.name) 

"""
```
mass(name::String or atom::Atom or Vector{Atom})
```

Returns the mass of an atom given its name, or `Atom` structure, or the total mass of a vector of `Atom`s. 

### Example

```julia-repl
julia> atoms = [ Atom(name="NT3"), Atom(name="CA") ];

julia> mass(atoms[1])
14.0067

julia> mass("CA")
12.011

julia> mass(atoms)
26.017699999999998

```

"""
mass(name::String) = elements[element_index(name)].mass
mass(atom::Atom) = mass(atom.name)
mass(atoms::AbstractVector{Atom}) = sum(mass(atoms[i]) for i in eachindex(atoms))

