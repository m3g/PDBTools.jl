#
# Function that reads atom information from PDB or mmCIF files
#
function read_atom(
    record::String;
    mmCIF::Bool=false,
    mmCIF_fields::Indexes_mmCIF_fields=Indexes_mmCIF_fields()
)
    atom = if mmCIF
        read_atom_mmCIF(record, mmCIF_fields)
    else
        read_atom_PDB(record)
    end
    return atom
end

function parse_number(::Type{T}, string, range) where {T<:AbstractFloat}
    parse(T, @view(string[range]))
end

function parse_number(::Type{T}, string, range) where {T<:Integer}
    s = @view(string[range])
    i = tryparse(Int, s)
    !isnothing(i) && return i
    # try to parse as hexadecimal number
    i = tryparse(Int, s, base=16)
    !isnothing(i) && return i
    error("Could not read integer from string: \"$s\"")
end

function parse_string(string, range)
    first_char = findfirst(>(' '), @view(string[range]))
    if isnothing(first_char)
        return @view(string[first(range):first(range)])
    end
    first_char = first(range) + first_char - 1
    last_char = first(range) + findlast(>(' '), @view(string[range])) - 1
    return @view(string[first_char:last_char])
end

# read atom from PDB file
function read_atom_PDB(record::String)
    N = length(record)
    if N < 6
        return nothing
    end
    if !(parse_string(record, 1:4) == "ATOM" || parse_string(record, 1:6) == "HETATM")
        return nothing
    end
    atom = Atom()
    atom.name = parse_string(record, 13:16)
    atom.resname = parse_string(record, 17:21)
    atom.chain = parse_string(record, 22:22)
    if isempty(atom.chain)
        atom.chain = "0"
    end
    atom.index = 1
    atom.index_pdb = parse_number(Int, record, 7:11)
    atom.resnum = parse_number(Int, record, 23:26)
    atom.x = parse_number(Float64, record, 31:38)
    atom.y = parse_number(Float64, record, 39:46)
    atom.z = parse_number(Float64, record, 47:54)
    atom.beta = parse_number(Float64, record, 61:66)
    atom.occup = parse_number(Float64, record, 56:60)
    atom.model = 1
    if N < 76
        atom.segname = "-"
    else
        atom.segname = parse_string(record, 73:76)
    end
    return atom
end

# read atom from mmCIF file
function read_atom_mmCIF(record::String, mmCIF_fields::Indexes_mmCIF_fields=Indexes_mmCIF_fields())
    if length(record) < 6 || !(record[1:4] == "ATOM" || record[1:6] == "HETATM")
        return nothing
    end
    atom = Atom()
    mmcif_data = split(record)
    atom.index = 1
    try
        atom.index_pdb = parse(Int, mmcif_data[mmCIF_fields.index])
    catch
        atom.index_pdb = 0
    end
    atom.name = mmcif_data[mmCIF_fields.name]
    atom.resname = mmcif_data[mmCIF_fields.resname]
    atom.chain = mmcif_data[mmCIF_fields.chain]
    try
        atom.resnum = parse(Int, mmcif_data[mmCIF_fields.resnum])
    catch
        atom.resnum = 0
    end
    try
        atom.segname = mmcif_data[mmCIF_fields.segname]
    catch
        atom.segname = ""
    end
    atom.x = parse(Float64, mmcif_data[mmCIF_fields.x])
    atom.y = parse(Float64, mmcif_data[mmCIF_fields.y])
    atom.z = parse(Float64, mmcif_data[mmCIF_fields.z])
    atom.beta = parse(Float64, mmcif_data[mmCIF_fields.beta])
    atom.occup = parse(Float64, mmcif_data[mmCIF_fields.occup])
    atom.model = 1
    return atom
end