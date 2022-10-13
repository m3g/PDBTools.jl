#
# Returns the index in the natural_aminoacids list of the amino acid name given
#
function which_natural_aminoacid(atom::Atom)

    # If the residue name doesn't have at least three letters, this is not a protein atom
    if length(atom.resname) < 3
        return 0
    end

    i = 0
    # First, check if there is a 4-letter code that matches the atom name
    i = findfirst(aa -> aa.three_letter_code == atom.resname, natural_aminoacids)
    if !isnothing(i)
        return i
    end

    # If not, try the same but removing the first char, to take into account
    # alternate conformations, such as "AGLY", "BGLY", etc.
    l = length(atom.resname)
    name = atom.resname[l-2:l]

    # Return the index of this amino acid in the amino acid list, or nothing
    i = findfirst(aa -> aa.three_letter_code == name, natural_aminoacids)
    if isnothing(i)
        return 0
    end

    return i
end
