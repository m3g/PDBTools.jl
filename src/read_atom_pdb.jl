#
# Function that reads atom information from PDB or mmCIF files
#
function read_atom(
    record::String;
    mmCIF::Bool=false,
    mmCIF_fields::indices_mmCIF_fields=indices_mmCIF_fields()
)
    atom = if mmCIF
        read_atom_mmCIF(record, mmCIF_fields)
    else
        read_atom_PDB(record)
    end
    return atom
end

function _parse(::Type{T}, string, range=nothing) where {T<:AbstractFloat}
    s = if isnothing(range)
        string
    else
        @view(string[range[begin]:min(range[end],length(string))])
    end
    x = tryparse(T, s)
    return isnothing(x) ? zero(T) : x
end

function _parse(::Type{T}, string, range=nothing; alt=nothing) where {T<:Integer}
    s = if isnothing(range)
        string
    else
        @view(string[range[begin]:min(range[end],length(string))])
    end
    i = tryparse(Int, s)
    !isnothing(i) && return i
    # try to parse as hexadecimal number
    i = tryparse(Int, s, base=16)
    !isnothing(i) && return i
    if isnothing(alt)
        error("Could not read integer from string: \"$s\"")
    else
        return alt
    end
end

function _parse(::Type{S}, string, range=nothing; alt=" ") where {S<:AbstractString}
    s, range = if isnothing(range)
        string, 1:length(string)
    else
        @view(string[range[begin]:min(range[end],length(string))]), range
    end
    length(range) > 0 || return alt
    first_char = findfirst(>(' '), s)
    isnothing(first_char) && return alt
    first_char = first(range) + first_char - 1
    last_char = first(range) + findlast(>(' '), s) - 1
    return @view(string[first_char:last_char])
end

# read atom from PDB file
function read_atom_PDB(record::String)
    if !startswith(record, r"ATOM|HETATM")
        return nothing
    end
    atom = Atom()
    atom.name = _parse(String, record, 13:16; alt="X")
    atom.resname = _parse(String, record, 17:21; alt="XXX")
    atom.chain = _parse(String, record, 22:22; alt="X")
    atom.index = 1
    atom.index_pdb = _parse(Int, record, 7:11; alt=-1)
    atom.resnum = _parse(Int, record, 23:26)
    atom.x = _parse(Float64, record, 31:38)
    atom.y = _parse(Float64, record, 39:46)
    atom.z = _parse(Float64, record, 47:54)
    atom.occup = _parse(Float64, record, 56:60)
    atom.beta = _parse(Float64, record, 61:66)
    atom.model = 1
    atom.segname = _parse(String, record, 73:76; alt="-")
    atom.pdb_element = _parse(String, record, 77:78, alt="X")
    atom.charge = _parse(String, record, 79:80, alt=nothing)
    return atom
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

