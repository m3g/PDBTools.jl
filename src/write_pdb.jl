"""
    write_pdb(atoms::Vector{Atom}, filename, selection; header=:auto, footer=:auto)

Write a PDB file with the atoms in `atoms` to `filename`. The `selection` argument is a string
that can be used to select a subset of the atoms in `atoms`. For example, `write_pdb(atoms, "test.pdb", "name CA")`.

The `header` and `footer` arguments can be used to add a header and footer to the PDB file. If `header` is `:auto`,
then a header will be added with the number of atoms in `atoms`. If `footer` is `:auto`, then a footer will be added
with the "END" keyword. Either can be set to `nothing` if no header or footer is desired.

"""
function write_pdb(atoms::AbstractVector{Atom}, filename::String, selection::String; header=:auto, footer=:auto)
    query = parse_query(selection)
    write_pdb(atoms, filename, only=atom -> apply_query(query, atom); header, footer)
end

function write_pdb(atoms::AbstractVector{Atom}, filename::String; only::Function=all, header=:auto, footer=:auto)
    file = open(expanduser(filename), "w")
    if header == :auto
        curr_date = Dates.format(Dates.today(), "dd-u-yy")
        header = "PDBTools.jl - $(length(atoms)) atoms"
        println(file, @sprintf "%-10s%-40s%9s" "HEADER" header curr_date)
    elseif header !== nothing
        println(file, header)
    end
    for atom in atoms
        if only(atom)
            println(file, write_pdb_atom(atom))
        end
    end
    if footer == :auto
        println(file, "END")
    elseif footer !== nothing
        println(file, footer)
    end
    close(file)
end

@testitem "write_pdb" begin
    using PDBTools
    using DelimitedFiles
    pdb = read_pdb(PDBTools.SMALLPDB)
    tmpfile = tempname()*".pdb"
    write_pdb(pdb, tmpfile)
    @test isfile(tmpfile)
    f1 = readdlm(PDBTools.SMALLPDB, '\n', header=true)
    f2 = readdlm(tmpfile, '\n', header=true)
    @test f1[1] == f2[1]
end
