"""
    readPDB(pdbfile::String, selection::String)
    readPDB(pdbfile::String; only::Function = all)

    readPDB(pdbdata::IOBuffer, selection::String)
    readPDB(pdbdata::IOBuffer; only::Function = all)

Reads a PDB file and stores the data in a vector of type `Atom`. 

If a selection is provided, only the atoms matching the selection will be read. 
For example, `resname ALA` will select all the atoms in the residue ALA.

If the `only` function keyword is provided, only the atoms for which `only(atom)` is true will be read.

### Examples

```julia-repl
julia> protein = readPDB("../test/structure.pdb")
   Array{Atoms,1} with 62026 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
                                                       ⋮ 
   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025
   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026

julia> ALA = readPDB("../test/structure.pdb","resname ALA")
   Array{Atoms,1} with 72 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
                                                       ⋮ 
    1339    C     ALA     A       95       95   14.815   -3.057   -5.633  0.00  1.00     1    PROT      1339
    1340    O     ALA     A       95       95   14.862   -2.204   -6.518  0.00  1.00     1    PROT      1340

julia> ALA = readPDB("../test/structure.pdb", only = atom -> atom.resname == "ALA")
   Array{Atoms,1} with 72 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
                                                       ⋮ 
    1339    C     ALA     A       95       95   14.815   -3.057   -5.633  0.00  1.00     1    PROT      1339
    1340    O     ALA     A       95       95   14.862   -2.204   -6.518  0.00  1.00     1    PROT      1340
```

"""
function readPDB end

function readPDB(file::Union{String,IOBuffer}, selection::String)
    query = parse_query(selection)
    return readPDB(file, only=atom -> apply_query(query, atom))
end

function readPDB(pdbdata::IOBuffer; only::Function=all)
    atoms = _parse_pdb(pdbdata, only)
    return atoms
end

function readPDB(file::String; only=all)
    atoms = open(expanduser(file), "r") do f
        _parse_pdb(f, only)
    end
    return atoms
end

function _parse_pdb(
    pdbdata::Union{IOStream, IOBuffer}, 
    only::Function, 
)
    imodel = 1
    atoms = Atom{Nothing}[]
    lastatom = Atom(; index=Int32(0), residue=Int32(0))
    for line in eachline(pdbdata)
        if occursin("END", line)
            imodel = imodel + 1
        end
        if startswith(line, r"ATOM|HETATM")
            atom = read_atom_PDB(line)
            atom.index = index(lastatom) + 1
            atom.model = imodel
            if !same_residue(atom, lastatom)
                atom.residue = residue(lastatom) + 1
            end
            if only(atom)
                push!(atoms, atom)
            end
            lastatom = atom
        end
    end
    seekstart(pdbdata)
    if length(atoms) == 0
        throw(ArgumentError("""\n 
            Could not find any atom in PDB file matching the selection. 

        """))
    end
    return atoms
end

# read atom from PDB file
function read_atom_PDB(record::String, atom = Atom(;index=Int32(1), model=Int32(1), custom=nothing))
    !startswith(record, r"ATOM|HETATM") && return nothing
    inds_and_names = (
        (1, Val(:name)), 
        (2, Val(:resname)), 
        (3, Val(:chain)), 
        (4, Val(:index_pdb)), 
        (5, Val(:resnum)), 
        (6, Val(:x)), 
        (7, Val(:y)), 
        (8, Val(:z)), 
        (9, Val(:occup)), 
        (10,Val(:beta)), 
        (11,Val(:segname)), 
        (12,Val(:pdb_element)), 
        (13,Val(:charge)),
    )
    @views field_values = (
        record[13:16], # :name
        record[17:21], # :resname
        record[22:22], # :chain
        record[7:11],  # :index_pdb
        record[23:26], # :resnum
        record[31:38], # :x
        record[39:46], # :y
        record[47:54], # :z
        length(record) >= 60 ? record[56:60] : record[56:end], # :occup
        length(record) >= 66 ? record[61:66] : record[61:end], # :beta
        length(record) >= 76 ? record[73:76] : record[73:end], # :segname
        length(record) >= 78 ? record[77:78] : record[77:end], # :pdb_element
        length(record) >= 80 ? record[79:80] : record[79:end], # :charge 
    )
    setfield_recursive!(atom, field_values, inds_and_names)
    return atom
end


@testitem "readPDB" begin
    pdb_file = "$(@__DIR__)/../test/structure.pdb"
    atoms = readPDB(pdb_file, "protein and name CA")
    @test length(atoms) == 104
    pdbdata = read(pdb_file, String)
    atoms = readPDB(IOBuffer(pdbdata), "protein and name CA")
    @test length(atoms) == 104
end

@testitem "read_atom_pdb" begin
    line = "HETATM    1  O1  BGL    1       0.665   1.214  -0.259  1.00  0.00"
    a = PDBTools.read_atom(line)
    @test a.index == 1
    @test a.name == "O1"
    @test a.resname == "BGL"
    @test a.chain == "X"
    @test a.resnum == 1
    @test a.residue == 0
    @test (a.x, a.y, a.z) == (0.665, 1.214, -0.259)
    @test (a.occup, a.beta) == (1.0, 0.0)
    @test a.model == 1
    @test a.segname == "-"
    @test a.index_pdb == 1
    @test a.pdb_element == "X"
    @test isnothing(a.charge)

    line = "HETATM    1  O1  AGL     1       0.803   1.186  -0.211  1.00  0.00            1+  "
    a = PDBTools.read_atom(line)
    @test a.index == 1
    @test a.name == "O1"
    @test a.resname == "AGL"
    @test a.chain == "X"
    @test a.resnum == 1
    @test a.residue == 0
    @test (a.x, a.y, a.z) == (0.803, 1.186, -0.211)
    @test (a.occup, a.beta) == (1.0, 0.0)
    @test a.model == 1
    @test a.segname == "-"
    @test a.index_pdb == 1
    @test a.pdb_element == "X"
    @test a.charge == "1+"

    line = "ATOM  *****  H2  GLYCD4301      11.014  36.823  40.115  1.00  0.00      GLYC H"
    a = PDBTools.read_atom(line)
    @test a.index == 1
    @test a.name == "H2"
    @test a.resname == "GLYC"
    @test a.chain == "D"
    @test a.resnum == 4301
    @test a.residue == 0
    @test (a.x, a.y, a.z) == (11.014, 36.823, 40.115)
    @test (a.occup, a.beta) == (1.0, 0.0)
    @test a.model == 1
    @test a.segname == "GLYC"
    @test a.index_pdb == -1
    @test a.pdb_element == "H"
    @test isnothing(a.charge)

    line = "ATOM  3ce2c  LP2 WLS Cffff     376.512 638.670  16.990  0.00  0.00         C"
    a = PDBTools.read_atom(line)
    @test a.index == 1
    @test a.index_pdb == 249388
    @test a.resnum == 65535
end


 