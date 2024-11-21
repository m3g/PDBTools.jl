"""
    getseq(AbstractVector{<:Atom} or filename; selection, code)

Returns the sequence of aminoacids from the vector of atoms or file name. Selections may be applied. Code defines if the output will be a one-letter, three-letter or full-residue name array.

### Example

```jldoctest
julia> using PDBTools

julia> protein = read_pdb(PDBTools.TESTPDB);

julia> getseq(protein,"residue < 3")
2-element Vector{String}:
 "A"
 "C"

julia> getseq(protein,"residue < 3", code=2)
2-element Vector{String}:
 "ALA"
 "CYS"

julia> getseq(protein,"residue < 3",code=3)
2-element Vector{String}:
 "Alanine"
 "Cysteine"

```

"""
function getseq(atoms::AbstractVector{<:Atom}, selection::String; code::Int=1)
    query = parse_query(selection)
    return getseq(atoms, only=atom -> apply_query(query, atom), code=code)
end

function getseq(atoms::AbstractVector{<:Atom}; only=isprotein, code::Int=1)
    seq = String[]
    for residue in eachresidue(atoms)
        # If any atom of this residue is in the selection, add it
        consider = false
        for at in residue
            if only(at)
                consider = true
                break
            end
        end
        if consider
            if isprotein(residue)
                code == 1 && push!(seq, oneletter(residue.name))
                code == 2 && push!(seq, threeletter(residue.name))
                code == 3 && push!(seq, residuename(residue.name))
            else
                code == 1 && push!(seq, name(residue)[1:1])
                code == 2 && push!(seq, name(residue))
                code == 3 && push!(seq, name(residue))
            end
        end
    end
    return seq
end

# From the file name
function getseq(file::String, selection::String; code::Int=1)
    atoms = read_pdb(file)
    return getseq(atoms, selection, code=code)
end

getseq(file::String; only=all, code::Int=1) = getseq(read_pdb(file), only=only, code=code)

@testitem "getseq" begin
    using PDBTools
    ats = read_pdb(PDBTools.TESTPDB)
    s = getseq(ats, "residue < 3")
    @test s == ["A", "C"]
    s = getseq(ats, "residue < 3"; code=2)
    @test s == ["ALA", "CYS"]
    s = getseq(ats, "residue < 3"; code=3)
    @test s == ["Alanine", "Cysteine"]
    wat = select(ats, "water")
    s = getseq(wat, "resnum < 3")
    @test s == ["T", "T", "T", "T"] 
    s = getseq(wat, "resnum < 3"; code=2)
    @test s == ["TIP3", "TIP3", "TIP3", "TIP3"] 
    s = getseq(wat, "resnum < 3"; code=3)    
    @test s == ["TIP3", "TIP3", "TIP3", "TIP3"] 
    s = getseq(PDBTools.TESTPDB, "residue < 3")
    s = getseq(ats, "residue < 3")
    @test s == ["A", "C"]
    s = getseq(PDBTools.TESTPDB; only = at -> resnum(at) == 1)
    @test s == ["A", "C", "S", "T", "T", "T"]
end
