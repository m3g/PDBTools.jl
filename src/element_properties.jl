#
# Try to deduce element (atomic number) from the atom name.
# Returns 0 instead.
#
function atomic_number(name :: String)

  # If the residue name doesn't have at least three letters, this is not a protein atom
  len = length(name)
  if len < 1
    return 0
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
    return elements[i].atomic_number
 
  # If not, check if the first character is a number, remove it and try again
  else
    try 
      parse(Int64,name[1:1])
      newname = name[2:length(name)]
      i = findfirst( el -> el.pdb_name == newname, elements )
      if i != nothing
        return elements[i].atomic_number
      else
        return 0
      end
    catch
      return 0
    end
  end

end

atomic_number(atom :: Union{Atom,MutableAtom}) = atomic_number(atom.name)

#
# Retrieve mass of atom
#
function mass(atom :: Union{Atom,MutableAtom})
  i = atomic_number(atom.name)
  elements[i].mass
end

#
# Retrieve name of atom
#
function name(atom :: Union{Atom,MutableAtom})
  i = atomic_number(atom.name)
  elements[i].name
end

#
# Retrieve name of atom
#
function code(atom :: Union{Atom,MutableAtom})
  i = atomic_number(atom.name)
  elements[i].code
end

