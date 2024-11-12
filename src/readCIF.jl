"""
    readCIF(mmCIF_file::String, selection::String)
    readCIF(mmCIF_file::String; only::Function = all)

    readCIF(mmCIF_data::IOBuffer, selection::String)
    readCIF(mmCIF_data::IOBuffer; only::Function = all)

Reads a mmCIF file and stores the data in a vector of type `Atom`. 

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

function readCIF(file::Union{String,IOBuffer}, selection::String; kargs...)
    query = parse_query(selection)
    return readCIF(file, only=atom -> apply_query(query, atom); kargs...)
end

function readCIF(cifdata::IOBuffer; only::Function=all, kargs...)
    atoms = _parse_mmCIF(cifdata; only, kargs...)
    return atoms
end

function readCIF(file::String; only=all, kargs...)
    atoms = open(expanduser(file), "r") do f
        _parse_mmCIF(f; only, kargs...)
    end
    return atoms
end

#_parse_mmCIF(args...;kargs...) = _parse_mmCIF_base(args...;kargs...)
#_parse_mmCIF(args...;kargs...) = _parse_mmCIF_csv(args...;kargs...)
_parse_mmCIF(args...;kargs...) = _parse_mmCIF_eachsplit(args...;kargs...)

function _parse_mmCIF_base(
    cifdata::Union{IOStream,IOBuffer};
    only::Function,
    memory_available::Real = 0.5,
    stop_at=nothing,
)
    _atom_site_field_inds = Dict{String, Int}()
    ifield = 0
    index = 0
    iresidue = 0
    atoms = LightAtom[]
    lastatom = LightAtom()
    for (iline, line) in enumerate(eachline(cifdata))
        # Reading the headers of the _atom_site loop
        if occursin("_atom_site.", line)
            field_end = findfirst(<=(' '), line) 
            if isnothing(field_end)
                field_end = length(line) + 1
            end
            field = @view(line[12:field_end-1])
            ifield += 1
            _atom_site_field_inds[field] = ifield
        end
        atom = read_atom_mmCIF_base(line, _atom_site_field_inds)
        if !isnothing(atom)
            index += 1
            atom.index = index
            if !same_residue(atom, lastatom)
                iresidue += 1
            end
            atom.residue = iresidue
            if only(atom)
                push!(atoms, atom)
            end
            lastatom = atom
        end
        if !isnothing(stop_at) && (length(atoms) >= stop_at)
            break
        end
        if mod(iline, 1000) == 0
            if Sys.free_memory() < (1 - memory_available) * Sys.total_memory()
                @warn """\n
                    Memory limit reached. $(length(atoms)) atoms read so far will be returned.
                    Size of the atoms array: $(round(Base.summarysize(atoms) / 1024^3; digits=3)) MB

                """ _file = nothing _line = nothing
                return atoms
            end
        end
    end
    seekstart(cifdata)
    if length(atoms) == 0
        throw(ArgumentError("""\n 
            Could not find any atom in mmCIF file matching the selection. 

        """))
    end
    return atoms
end
#=
data_all
loop_
_atom_site.group_PDB             # ATOM 
_atom_site.id                    # 1 index_pdb
_atom_site.type_symbol           # N symbol?  
_atom_site.label_atom_id         # N name
_atom_site.label_alt_id          # . 
_atom_site.label_comp_id         # VAL resname
_atom_site.label_asym_id         # A chain
_atom_site.label_entity_id       # 1 ?
_atom_site.label_seq_id          # 1 resnum
_atom_site.pdbx_PDB_ins_code     # ? ? 
_atom_site.Cartn_x               # 6.204 x
_atom_site.Cartn_y               # 16.869 y
_atom_site.Cartn_z               # 4.854 z
_atom_site.occupancy             # 1.00 occup
_atom_site.B_iso_or_equiv        # 49.05 beta
_atom_site.pdbx_formal_charge    # ? charge
_atom_site.auth_seq_id           # 1
_atom_site.auth_comp_id          # VAL
_atom_site.auth_asym_id          # A
_atom_site.auth_atom_id          # N
_atom_site.pdbx_PDB_model_num    # 1 model
ATOM   1    N  N   . VAL A 1 1   ? 6.204   16.869  4.854   1.00 49.05 ? 1   VAL A N   1
ATOM   2    C  CA  . VAL A 1 1   ? 6.913   17.759  4.607   1.00 43.14 ? 1   VAL A CA  1
ATOM   3    C  C   . VAL A 1 1   ? 8.504   17.378  4.797   1.00 24.80 ? 1   VAL A C   1
ATOM   4    O  O   . VAL A 1 1   ? 8.805   17.011  5.943   1.00 37.68 ? 1   VAL A O   1
ATOM   5    C  CB  . VAL A 1 1   ? 6.369   19.044  5.810   1.00 72.12 ? 1   VAL A CB  1
ATOM   6    C  CG1 . VAL A 1 1   ? 7.009   20.127  5.418   1.00 61.79 ? 1   VAL A CG1 1
ATOM   7    C  CG2 . VAL A 1 1   ? 5.246   18.533  5.681   1.00 80.12 ? 1   VAL A CG2 1 
=#

# read atom from mmCIF file
function read_atom_mmCIF_base(record::String, _atom_site_field_inds)
    if !startswith(record, r"ATOM|HETATM")
        return nothing
    end
    atom = LightAtom(;index=1, segname="")
    mmcif_data = split(record)
    atom.index_pdb = _parse(Int32, mmcif_data[_atom_site_field_inds["id"]])
    atom.name = mmcif_data[_atom_site_field_inds["label_atom_id"]]
    atom.resname = mmcif_data[_atom_site_field_inds["label_comp_id"]]
    atom.chain = mmcif_data[_atom_site_field_inds["label_asym_id"]]
    atom.resnum = _parse(Int32, mmcif_data[_atom_site_field_inds["label_seq_id"]])
    atom.x = _parse(Float32, mmcif_data[_atom_site_field_inds["Cartn_x"]])
    atom.y = _parse(Float32, mmcif_data[_atom_site_field_inds["Cartn_y"]])
    atom.z = _parse(Float32, mmcif_data[_atom_site_field_inds["Cartn_z"]])
    if haskey(_atom_site_field_inds, "B_iso_or_equiv")
        atom.beta = _parse(Float32, mmcif_data[_atom_site_field_inds["B_iso_or_equiv"]])
    end
    if haskey(_atom_site_field_inds, "occupancy")
        atom.occup = _parse(Float32, mmcif_data[_atom_site_field_inds["occupancy"]])
    end
    if haskey(_atom_site_field_inds, "pdbx_PDB_model_num")
        atom.model = _parse(Int32, mmcif_data[_atom_site_field_inds["pdbx_PDB_model_num"]])
    end
    if haskey(_atom_site_field_inds, "pdbx_formal_charge")
        atom.charge = _parse(Float32, mmcif_data[_atom_site_field_inds["pdbx_formal_charge"]])
    end
    return atom
end

#
# test with eachsplit
#
function _parse_mmCIF_eachsplit(
    cifdata::Union{IOStream,IOBuffer};
    only::Function,
    memory_available::Real=0.5,
    stop_at=nothing,
)
    _atom_symbol_for_cif_field = Dict{String, Tuple{DataType,Symbol}}(
        "id" => (Int32, :index_pdb),
        "label_atom_id" => (String7, :name),
        "label_comp_id" => (String7, :resname),
        "label_asym_id" => (String3, :chain),
        "label_seq_id" => (Int32, :resnum),
        "Cartn_x" => (Float32,:x),
        "Cartn_y" => (Float32,:y),
        "Cartn_z" => (Float32,:z),
        "occupancy" => (Float32,:occup),
        "B_iso_or_equiv" => (Float32,:beta),
        "pdbx_formal_charge" => (Float32,:charge),
        "pdbx_PDB_model_num" => (Int32,:model),
    )
    _atom_site_field_inds = Dict{String,Int}()
    ifield = 0
    atoms = LightAtom[]
    lastatom = LightAtom()
    _atom_field_columns = Vector{Tuple{Int,Tuple{DataType,Symbol}}}()
    local NCOLS, col_indices, col_types, col_field
    for line in eachline(cifdata)
        # Reading the headers of the _atom_site loop
        if occursin("_atom_site.", line)
            field_end = findfirst(<=(' '), line) 
            if isnothing(field_end)
                field_end = length(line) + 1
            end
            field = @view(line[12:field_end-1])
            ifield += 1
            _atom_site_field_inds[field] = ifield
        end
        # Header ended
        if startswith(line, r"ATOM|HETATM")
            for key in keys(_atom_symbol_for_cif_field)
                if haskey(_atom_site_field_inds, key)
                    push!(_atom_field_columns, (_atom_site_field_inds[key], _atom_symbol_for_cif_field[key]))
                end
            end
            sort!( _atom_field_columns; by = first)
            col_indices = NTuple{length(_atom_field_columns),Int}(first(el) for el in _atom_field_columns)
            col_types = NTuple{length(_atom_field_columns),DataType}(first(last(el)) for el in _atom_field_columns)
            col_field = NTuple{length(_atom_field_columns),Symbol}(last(last(el)) for el in _atom_field_columns)
            NCOLS = length(keys(_atom_site_field_inds))
            break
        end
    end
    seekstart(cifdata)
    for line in eachline(cifdata)
        !startswith(line, r"ATOM|HETATM") && continue
        # Assign the atom fields to each column
        atom = read_atom_mmCIF_eachsplit(
            Val(NCOLS), line, col_indices, col_types, col_field;
            lastatom,
        )
        lastatom = atom
        if only(atom)
            push!(atoms, atom)
        end
        if !isnothing(stop_at) && (length(atoms) >= stop_at)
            break
        end
        if mod(length(atoms), 1000) == 0
            if Sys.free_memory() < (1 - memory_available) * Sys.total_memory()
                @warn """\n
                    Memory limit reached. $(length(atoms)) atoms read so far will be returned.
                    Size of the atoms array: $(round(Base.summarysize(atoms) / 1024^3; digits=3)) MB

                """ _file = nothing _line = nothing
                return atoms
            end
        end
    end
    seekstart(cifdata)
    if length(atoms) == 0
        throw(ArgumentError("""\n 
            Could not find any atom in mmCIF file matching the selection. 

        """))
    end
    return atoms
end

function read_atom_mmCIF_eachsplit(
    ::Val{NCOLS}, record::String, col_indices, col_types, col_field;
    lastatom::AbstractAtom,
) where {NCOLS}
    fields = NTuple{NCOLS}(eachsplit(record))
    atom = LightAtom(; index = index(lastatom) + 1, residue = residue(lastatom))
    for (ifield, field) in enumerate(col_field)
        if col_types[ifield] == Float32
            setfield!(atom, field, parse(Float32, fields[col_indices[ifield]]))
        end
        if col_types[ifield] == Int32
            setfield!(atom, field, parse(Int32, fields[col_indices[ifield]]))
        end
        if col_types[ifield] == String7
            setfield!(atom, field, String7(fields[col_indices[ifield]]))
        end
        if col_types[ifield] == String3
            setfield!(atom, field, String3(fields[col_indices[ifield]]))
        end
    end
    #for (ifield, field) in enumerate(col_field)
    #    setfield!(atom, field, _parse(col_types[ifield], fields[col_indices[ifield]]))
    #end
    if !same_residue(atom, lastatom)
        atom.residue = residue(lastatom) + 1
    end
    return atom
end

#
# test with CSV.Rows
#
function _parse_mmCIF_csv(
    cifdata::Union{IOStream,IOBuffer};
    only::Function,
    memory_available::Real=0.5,
    stop_at=nothing,
    float_type::DataType=Float64,
)
    _atom_symbol_for_cif_field = Dict{String, Tuple{DataType,Symbol}}(
        "id" => (Int, :index_pdb),
        "label_atom_id" => (String, :name),
        "label_comp_id" => (String, :resname),
        "label_asym_id" => (String, :chain),
        "label_seq_id" => (Int, :resnum),
        "Cartn_x" => (float_type,:x),
        "Cartn_y" => (float_type,:y),
        "Cartn_z" => (float_type,:z),
        "occupancy" => (float_type,:occup),
        "B_iso_or_equiv" => (float_type,:beta),
        "pdbx_formal_charge" => (String,:charge),
        "pdbx_PDB_model_num" => (Int,:model),
    )
    _atom_site_field_inds = Dict{String,Int}()
    ifield = 0
    _atom_field_columns = Vector{Tuple{Int,Tuple{DataType,Symbol}}}()
    local col_types, col_field
    for line in eachline(cifdata)
        # Reading the headers of the _atom_site loop
        if occursin("_atom_site.", line)
            field_end = findfirst(<=(' '), line) 
            if isnothing(field_end)
                field_end = length(line) + 1
            end
            field = @view(line[12:field_end-1])
            ifield += 1
            _atom_site_field_inds[field] = ifield
        end
        if startswith(line, r"ATOM|HETATM")
            for key in keys(_atom_site_field_inds)
                if haskey(_atom_symbol_for_cif_field, key)
                    push!(_atom_field_columns, (_atom_site_field_inds[key], _atom_symbol_for_cif_field[key]))
                else
                    push!(_atom_field_columns, (_atom_site_field_inds[key], (String, :ignore)))
                end
            end
            sort!( _atom_field_columns; by = first)
            # this is used to identify atom lines later
            _atom_field_columns[1] = (1, (String, :atom))
            col_types = DataType[first(last(el)) for el in _atom_field_columns]
            col_field = Symbol[last(last(el)) for el in _atom_field_columns]
            break
        end
    end
    col_types = SVector{length(col_types),DataType}(col_types)
    seekstart(cifdata)
    atoms = read_atoms_mmCIF_csv(float_type,cifdata, col_types, col_field, only, memory_available, stop_at)
    if length(atoms) == 0
        throw(ArgumentError("""\n 
            Could not find any atom in mmCIF file matching the selection. 

        """))
    end
    return atoms
end

function read_atoms_mmCIF_csv(::Type{float_type}, cifdata, col_types, col_field, only, memory_available, stop_at) where {float_type}
    index = 0
    iresidue = 0
    lastatom = Atom{float_type}()
    atoms = Atom{float_type}[]
    for row in CSV.Rows(cifdata; 
        delim=' ', 
        stripwhitespace=true, 
        ignorerepeated=true, 
        header=col_field, 
        types=(i,name) -> col_types[i], 
        silencewarnings=true,
        reusebuffer=true,
    )
        if row.atom in ("ATOM","HETATM")
            atom = Atom{float_type}()
            index += 1
            for field in col_field
                if field in (:atom, :ignore)
                    continue
                end
                setfield!(atom, field, getproperty(row, field))
            end
            atom.index = index
            atom.segname = ""
            if !same_residue(atom, lastatom)
                iresidue += 1
            end
            atom.residue = iresidue
            if only(atom)
                push!(atoms, atom)
            end
            lastatom = atom
        end
        if !isnothing(stop_at) && (length(atoms) >= stop_at)
            break
        end
        if mod(length(atoms), 1000) == 0
            if Sys.free_memory() < (1 - memory_available) * Sys.total_memory()
                 @warn """\n
                     Memory limit reached. $(length(atoms)) atoms read so far will be returned.
                     Size of the atoms array: $(round(Base.summarysize(atoms) / 1024^3; digits=3)) MB

                 """ _file = nothing _line = nothing
                 return atoms
             end
         end
    end
    return atoms
end
