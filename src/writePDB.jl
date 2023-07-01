"""
    writePDB(atoms::Vector{Atom}, filename, selection)

Function that writes a PDB file given the vector of atoms, with
optional definition of a selection to be print.

"""
function writePDB(atoms::AbstractVector{Atom}, filename, selection)
    query = parse_query(selection)
    writePDB(atoms, filename, only=atom -> apply_query(query, atom))
end

function writePDB(atoms::AbstractVector{Atom}, filename; only=all)
    file = open(filename, "w")
    curr_date = Dates.format(Dates.today(), "dd-u-yy")
    header = "PDBTools.jl - $(length(atoms)) atoms"
    println(file, @sprintf "%-10s%-40s%9s" "HEADER" header curr_date)
    for atom in atoms
        if only(atom)
            println(file, write_atom(atom))
        end
    end
    println(file, "END")
    close(file)
end
