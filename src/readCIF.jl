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
function readCIF end

function readCIF(file::Union{String,IOBuffer}, selection::String)
    query = parse_query(selection)
    return readCIF(file, only=atom -> apply_query(query, atom))
end

function readCIF(cifdata::IOBuffer; only::Function=all)
    atoms = _parse_mmCIF(cifdata, only)
    return atoms
end

function readCIF(file::String; only=all)
    atoms = open(expanduser(file), "r") do f
        _parse_mmCIF(f, only)
    end
    return atoms
end

function _parse_mmCIF(cifdata::Union{IOStream,IOBuffer})
    _atom_site_field_inds = Dict{String, Int}(
        "label_atom_id" => 0,        # N name
        "label_comp_id" => 0,        # VAL resname
        "label_asym_id" => 0,        # A chain
        "label_entity_id" => 0,      # 1 resnum
        "Cartn_x" => 0,              # 6.204 x
        "Cartn_y" => 0,              # 16.869 y
        "Cartn_z" => 0,              # 4.854 z
        "occupancy" => 0,            # 1.00 occup
        "B_iso_or_equiv" => 0,       # 49.05 beta
        "pdbx_formal_charge" => 0,   # ? charge
        "pdbx_PDB_model_num" => 0,   # 1 model
    )
    ifield = 0
    atoms = Atom[]
    for line in eachline(cifdata)
        # Reading the headers of the _atom_site loop
        if occursin("_atom_site.", line)
            field = @view(line[12:findfirst(<=(' '), line)-1])
            if haskey(_atom_site_field_inds, field)
                ifield += 1
                _atom_site_field_inds[field] = ifield
            end
        end
        # Parsing ATOM or HETATM lines
        if occursin("ATOM", line) || occursin("HETATM", line)
            atom = read_atom_mmCIF(line, _atom_site_field_inds)
            if !isnothing(atom)
                push!(atoms, atom)
            end
        end
    end
    return 
end

#=
function _parse_cif(
    pdbdata::Union{IOStream, IOBuffer}, 
    only::Function, 
    mmCIF::Bool,
    mmCIF_fields::indices_mmCIF_fields
)
    natoms = 0
    index = 0
    imodel = 1
    iresidue = 1
    atoms = Atom[]
    local lastatom
    for line in eachline(pdbdata)
        if occursin("END", line)
            imodel = imodel + 1
        end
        atom = read_atom(line, mmCIF=mmCIF, mmCIF_fields=mmCIF_fields)
        if !isnothing(atom)
            index = index + 1
            atom.index = index
            atom.model = imodel
            if index > 1
                if !same_residue(atom, lastatom)
                    iresidue += 1
                end
            end
            atom.residue = iresidue
            if only(atom)
                natoms = natoms + 1
                push!(atoms, atom)
            end
            lastatom = atom
        end
    end
    seekstart(pdbdata)
    if natoms == 0
        error(" Could not find any atom in PDB file matching the selection. ")
    end
    return atoms
end

function Base.show(io::IO, atoms::AbstractVector{Atom})
    println(io, " Structure file with ", length(atoms), " atoms. ")
end

@testitem "readPDB" begin
    pdb_file = "$(@__DIR__)/../test/structure.pdb"
    atoms = readPDB(pdb_file, "protein and name CA")
    @test length(atoms) == 104
    pdbdata = read(pdb_file, String)
    atoms = readPDB(IOBuffer(pdbdata), "protein and name CA")
    @test length(atoms) == 104
end
=#

