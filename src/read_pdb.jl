"""
    read_pdb(pdbfile::String, selection::String)
    read_pdb(pdbfile::String; only::Function = all)

    read_pdb(pdbdata::IOBuffer, selection::String)
    read_pdb(pdbdata::IOBuffer; only::Function = all)

Reads a PDB file and stores the data in a vector of type `Atom`. 

If a selection is provided, only the atoms matching the selection will be read. 
For example, `resname ALA` will select all the atoms in the residue ALA.

If the `only` function keyword is provided, only the atoms for which `only(atom)` is true will be read.

### Examples

```jldoctest
julia> using PDBTools

julia> protein = read_pdb(PDBTools.TESTPDB)
   Vector{Atom{Nothing}} with 62026 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
                                                       ⋮
   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  1.00  0.00     1    WAT2     62024
   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  1.00  0.00     1    WAT2     62025
   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  1.00  0.00     1    WAT2     62026

julia> ALA = read_pdb(PDBTools.TESTPDB,"resname ALA")
   Vector{Atom{Nothing}} with 72 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
                                                       ⋮
    1338  HB3     ALA     A       95       95   11.464   -3.278   -4.953  1.00  0.00     1    PROT      1338
    1339    C     ALA     A       95       95   14.815   -3.057   -5.633  1.00  0.00     1    PROT      1339
    1340    O     ALA     A       95       95   14.862   -2.204   -6.518  1.00  0.00     1    PROT      1340

julia> ALA = read_pdb(PDBTools.TESTPDB, only = atom -> atom.resname == "ALA")
   Vector{Atom{Nothing}} with 72 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
                                                       ⋮
    1338  HB3     ALA     A       95       95   11.464   -3.278   -4.953  1.00  0.00     1    PROT      1338
    1339    C     ALA     A       95       95   14.815   -3.057   -5.633  1.00  0.00     1    PROT      1339
    1340    O     ALA     A       95       95   14.862   -2.204   -6.518  1.00  0.00     1    PROT      1340
```

"""
function read_pdb end

function read_pdb(file::Union{String,IOBuffer}, selection::String)
    query = parse_query(selection)
    return read_pdb(file, only=atom -> apply_query(query, atom))
end

function read_pdb(pdbdata::IOBuffer; only::Function=all)
    atoms = _parse_pdb(pdbdata, only)
    return atoms
end

function read_pdb(file::String; only=all)
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
    lastatom = Atom{Nothing}()
    for line in eachline(pdbdata)
        if occursin("END", line)
            imodel = imodel + 1
        end
        if startswith(line, r"ATOM|HETATM")
            atom = read_atom_pdb(line, lastatom, imodel)
            atom.model = imodel
            only(atom) && push!(atoms, atom)
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
function read_atom_pdb(record::String, lastatom::Atom=Atom{Nothing}(), imodel::Int=1)
    atom = Atom{Nothing}(;index=index(lastatom) + 1, residue=residue(lastatom), model=imodel)
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
        length(record) >= 80 ? _parse_charge(record[79:80]) : _parse_charge(record[79:end]), # :charge 
    )
    _fast_setfield!(atom, field_values, inds_and_names)
    if !same_residue(atom, lastatom)
        atom.residue = residue(lastatom) + 1
    end
    return atom
end

@testitem "read_pdb" begin
    pdb_file = PDBTools.TESTPDB
    atoms = read_pdb(pdb_file, "protein and name CA")
    @test length(atoms) == 104
    pdbdata = read(pdb_file, String)
    atoms = read_pdb(IOBuffer(pdbdata), "protein and name CA")
    @test length(atoms) == 104
end

@testitem "read_atom_pdb" begin
    line = "HETATM    1  O1  BGL    1       0.665   1.214  -0.259  1.00  0.00"
    a = PDBTools.read_atom_pdb(line)
    @test a.index == 1
    @test a.name == "O1"
    @test a.resname == "BGL"
    @test a.chain == "X"
    @test a.resnum == 1
    @test a.residue == 1
    @test (a.x, a.y, a.z) == (0.665f0, 1.214f0, -0.259f0)
    @test (a.occup, a.beta) == (1.0, 0.0)
    @test a.model == 1
    @test a.segname == "X"
    @test a.index_pdb == 1
    @test a.pdb_element == "X"
    @test a.charge == 0.0f0

    line = "HETATM    1  O1  AGL     1       0.803   1.186  -0.211  1.00  0.00            1+  "
    a = PDBTools.read_atom_pdb(line)
    @test a.index == 1
    @test a.name == "O1"
    @test a.resname == "AGL"
    @test a.chain == "X"
    @test a.resnum == 1
    @test a.residue == 1
    @test (a.x, a.y, a.z) == (0.803f0, 1.186f0, -0.211f0)
    @test (a.occup, a.beta) == (1.0, 0.0)
    @test a.model == 1
    @test a.segname == "X"
    @test a.index_pdb == 1
    @test a.pdb_element == "X"
    @test a.charge == 1.0f0

    line = "HETATM    1  O1  AGL     1       0.803   1.186  -0.211  1.00  0.00            +1  "
    a = PDBTools.read_atom_pdb(line)
    @test a.charge == 1.0f0

    line = "HETATM    1  O1  AGL     1       0.803   1.186  -0.211  1.00  0.00            -1  "
    a = PDBTools.read_atom_pdb(line)
    @test a.charge == -1.0f0

    line = "HETATM    1  O1  AGL     1       0.803   1.186  -0.211  1.00  0.00            1-  "
    a = PDBTools.read_atom_pdb(line)
    @test a.charge == -1.0f0

    line = "ATOM  *****  H2  GLYCD4301      11.014  36.823  40.115  1.00  0.00      GLYC H"
    a = PDBTools.read_atom_pdb(line)
    @test a.index == 1
    @test a.name == "H2"
    @test a.resname == "GLYC"
    @test a.chain == "D"
    @test a.resnum == 4301
    @test a.residue == 1
    @test (a.x, a.y, a.z) == (11.014f0, 36.823f0, 40.115f0)
    @test (a.occup, a.beta) == (1.0, 0.0)
    @test a.model == 1
    @test a.segname == "GLYC"
    @test a.index_pdb == 0
    @test a.pdb_element == "H"
    @test a.charge == 0

    line = "ATOM  3ce2c  LP2 WLS Cffff     376.512 638.670  16.990  0.00  0.00         C"
    a = PDBTools.read_atom_pdb(line)
    @test a.index == 1
    @test a.index_pdb == 249388
    @test a.resnum == 65535
end
