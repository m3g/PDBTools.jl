#
# Retrive index of element in elements list from name. Returns 1 (element "X" of list) of not found
#
function element_index(name::String)

  # If the residue name doesn't have at least three letters, this is not a protein atom
  len = length(name)
  if len < 1
    return 1
  end

  # Return the index of this amino acid in the element, or nothing
  i = nothing
  l = len
  while i == nothing && l >= 0  
    i = findfirst( el -> el.pdb_name == name[1:l], elements )
    l = l - 1 
  end

  # If found, return the atomic number
  if i != nothing
    return i
 
  # If not, check if the first character is a number, remove it and try again
  else
    try 
      parse(Int, name[1:1])
      newname = name[2:length(name)]
      i = findfirst( el -> el.pdb_name == newname, elements )
      if i != nothing
        return i
      else
        return 1
      end
    catch
      return 1
    end
  end

end

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

