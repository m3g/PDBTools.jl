"""
    edit!(atoms::Vector{Atom})

Opens a temporary PDB file in which the fields of the vector of atoms can be edited.   

"""
function edit!(atoms::AbstractVector{Atom})
    tmp_file_name = tempname()
    write_pdb(tmp_file_name, atoms)
    InteractiveUtils.edit(tmp_file_name)
    read_again = read_pdb(tmp_file_name)
    resize!(atoms, length(read_again))
    @. atoms = read_again
    rm(tmp_file_name)
end
