"""
    readPDB(filename; selection)

Reads a PDB file and stores the data in a vector of type `Atom`.

### Example

```julia-repl
julia> protein = readPDB("../test/structure.pdb")
   Array{Atoms,1} with 62026 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
                                                       ⋮ 
   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  0.00  1.00     1    WAT2     62024
   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025
   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026

julia> protein = readPDB("../test/structure.pdb","resname = ALA")
   Array{Atoms,1} with 72 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
                                                       ⋮ 
    1338  HB3     ALA     A       95       95   11.464   -3.278   -4.953  0.00  1.00     1    PROT      1338
    1339    C     ALA     A       95       95   14.815   -3.057   -5.633  0.00  1.00     1    PROT      1339
    1340    O     ALA     A       95       95   14.862   -2.204   -6.518  0.00  1.00     1    PROT      1340
```

"""
function readPDB(file::String, selection::String)
    query = parse_query(selection)
    return readPDB(file, only = atom -> apply_query(query, atom))
end

function readPDB(file::String; only = all)

    # Check if structure is in mmCIF format
    mmCIF, mmCIF_fields = check_mmCIF(file)

    # Read file
    pdbfile = open(file, "r")
    natoms = 0
    index = 0
    imodel = 1
    iresidue = 1
    atoms = Atom[]
    local lastatom
    for line in eachline(pdbfile)
        if occursin("END", line)
            imodel = imodel + 1
        end
        atom = read_atom(line, mmCIF = mmCIF, mmCIF_fields = mmCIF_fields)
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
    close(pdbfile)
    if natoms == 0
        error(" Could not find any atom in PDB file matching the selection. ")
    end

    return atoms
end

function Base.show(io::IO, atoms::AbstractVector{Atom})
    println(io, " Structure file with ", length(atoms), " atoms. ")
end
