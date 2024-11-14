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

#
# Function that returns an ATOM line in PDB format
#
# requires Printf

function _align_name(name)
    name = strip(name)
    length(name) == 1 && return " $(name)  "
    length(name) == 2 && return " $(name) "
    length(name) == 3 && return " $(name)"
    return name
end
function _align_resname(resname)
    resname = strip(resname)
    length(resname) == 1 && return "  $(resname) "
    length(resname) == 2 && return " $(resname)  "
    length(resname) == 3 && return " $(resname) "
    return resname
end

function write_pdb_atom(atom::Atom)

    #ATOM      2  CA  GLY A   1      -1.774   6.778  32.054  1.00  0.08           C
    #COLUMNS        DATA  TYPE    FIELD        DEFINITION
    #-------------------------------------------------------------------------------------
    # 1 -  6        Record name   "ATOM  "
    # 7 - 11        Integer       serial       Atom  serial number.
    #13 - 16        Atom          name         Atom name.
    #17             Character     altLoc       Alternate location indicator.
    #18 - 20        Residue name  resName      Residue name.
    #22             Character     chainID      Chain identifier.
    #23 - 26        Integer       resSeq       Residue sequence number.
    #27             AChar         iCode        Code for insertion of residues.
    #31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    #39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    #47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    #55 - 60        Real(6.2)     occupancy    Occupancy.
    #61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    #73 - 76        String        segname      Segment identifier, left-justified (not default of PDB)
    #77 - 78        LString(2)    element      Element symbol string, right-justified.
    #79 - 80        LString(2)    charge       Charge  on the atom.

    name = strip(atom.name)
    l = length(name)
    if l == 1
        name = "  $(atom.name) "
    elseif l == 2
        name = " $(atom.name) "
    elseif l == 3
        name = " $(atom.name)"
    else
        name = atom.name
    end

    if abs(atom.x) > 999.999 || abs(atom.y) > 999.999 || abs(atom.y) > 999.999
        println(
            "Warning: coordinates of atom $(atom.index) " *
            "($(atom.name) $(atom.resname)$(atom.resnum) $(atom.chain)) " *
            "do not fit in PDB fixed format.",
        )
    end

    atom_data = (
            "ATOM",
            mod(atom.index,1048576),
            " ",
            _align_name(name),
            _align_resname(atom.resname),
            atom.chain,
            mod(atom.resnum, 65535), 
            "    ",
            atom.x,
            atom.y,
            atom.z,
            atom.occup,
            atom.beta,
            "      ",
            segname(atom),
            element(atom) === nothing ? "  " : element_symbol_string(atom),
    )

    line = if atom.index <= 99999 && atom.resnum <= 9999 # standard integer printing
        @sprintf("%-6s%5i%1s%4s%4s%1s%4i%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    elseif atom.index > 99999 && atom.resnum <= 9999 # Prints index in hexadecimal
        @sprintf("%-6s%5x%1s%-4s%4s%1s%4i%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    elseif atom.index <= 99999 && atom.resnum > 9999 # Prints resnum in hexadecimal code
        @sprintf("%-6s%5i%1s%-4s%4s%1s%4x%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    elseif atom.index > 99999 && atom.resnum > 9999 # Both hexadecimal
        @sprintf("%-6s%5x%1s%-4s%4s%1s%4x%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    end

    return line
end

@testitem "write_pdb_atom" begin
    using PDBTools
    pdb = read_pdb(PDBTools.SMALLPDB)
    @test PDBTools.write_pdb_atom(pdb[1]) == "ATOM      1  N   ALA A   1      -9.229 -14.861  -5.481  0.00  0.00      PROT N" 
    @test PDBTools.write_pdb_atom(pdb[2]) == "ATOM      2 1HT1 ALA A   1     -10.048 -15.427  -5.569  0.00  0.00      PROT H"
    pdb[1].index = 1000000
    @test PDBTools.write_pdb_atom(pdb[1]) == "ATOM  f4240  N   ALA A   1      -9.229 -14.861  -5.481  0.00  0.00      PROT N" 
    pdb[1].index = 1
    pdb[1].resnum = 1000000
    @test PDBTools.write_pdb_atom(pdb[1]) == "ATOM      1  N   ALA A424f      -9.229 -14.861  -5.481  0.00  0.00      PROT N" 
    pdb[1].index = 1000000
    @test PDBTools.write_pdb_atom(pdb[1]) == "ATOM  f4240  N   ALA A424f      -9.229 -14.861  -5.481  0.00  0.00      PROT N" 
end
