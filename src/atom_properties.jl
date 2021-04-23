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
      parse(Int64,name[1:1])
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

atomic_number(name::String) = elements[element_index(name)].atomic_number
atomic_number(atom::Atom) = atomic_number(atom.name)

element(name::String) = elements[element_index(name)].element
element(atom::Atom) = element(atom.name)

mass(name::String) = elements[element_index(name)].mass
mass(atom::Atom) = mass(atom.name)
mass(atoms::Vector{Atom}) = sum(at -> mass(at), atoms)

name(name::String) = elements[element_index(name)].name
name(atom::Atom) = name(atom.name)



