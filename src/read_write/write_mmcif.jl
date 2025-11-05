#=
loop_
_atom_site.group_PDB                  #  1 ATOM 
_atom_site.id                         #  2 index
_atom_site.type_symbol                #  3 element
_atom_site.label_atom_id              #  4 redundant: atom name (site.auth_atom_id)
_atom_site.label_alt_id               #  5 nothing
_atom_site.label_comp_id              #  6 redundant: residue name (site.auth_comp_id)
_atom_site.label_asym_id              #  7 redundant: chain id (site.auth_asym_id)
_atom_site.label_entity_id            #  8 entity_id: print 1 for protein, 2 for water, 3 for all others
_atom_site.label_seq_id               #  9 redundant: residue number (site.auth_seq_id)
_atom_site.pdbx_PDB_ins_code          # 10 nothing - print ?
_atom_site.Cartn_x                    # 11 x
_atom_site.Cartn_y                    # 12 y
_atom_site.Cartn_z                    # 13 z
_atom_site.occupancy                  # 14 occupancy
_atom_site.B_iso_or_equiv             # 15 bfactor
_atom_site.pdbx_formal_charge         # 16 charge: print ? if unknown
_atom_site.auth_seq_id                # 17 residue number   
_atom_site.auth_comp_id               # 18 residue name
_atom_site.auth_asym_id               # 19 chain id
_atom_site.auth_atom_id               # 20 atom name
_atom_site.pdbx_PDB_model_num         # 21 model number
ATOM   1    N  N   . VAL A 1 1   ? 6.204   16.869  4.854   1.00 49.05 ? 1   VAL A N   1
ATOM   2    C  CA  . VAL A 1 1   ? 6.913   17.759  4.607   1.00 43.14 ? 1   VAL A CA  1
ATOM   3    C  C   . VAL A 1 1   ? 8.504   17.378  4.797   1.00 24.80 ? 1   VAL A C   1
ATOM   4    O  O   . VAL A 1 1   ? 8.805   17.011  5.943   1.00 37.68 ? 1   VAL A O   1
ATOM   5    C  CB  . VAL A 1 1   ? 6.369   19.044  5.810   1.00 72.12 ? 1   VAL A CB  1
ATOM   6    C  CG1 . VAL A 1 1   ? 7.009   20.127  5.418   1.00 61.79 ? 1   VAL A CG1 1
ATOM   7    C  CG2 . VAL A 1 1   ? 5.246   18.533  5.681   1.00 80.12 ? 1   VAL A CG2 1 
=#
function _supported_write_cif_fields(field_assignment)
    # We need this to be indexable (and ordered dict) to keep the order of the fields
    # when writing the mmCIF file
    _atom_symbol_for_cif_field = OrderedDict{String,Tuple{DataType,Symbol}}(
        "id" => (Int32, :index_pdb), # Standard mmCIF
        "type_symbol" => (String7, :pdb_element), # Standard mmCIF
        "label_atom_id" => (String7, :name), # redundant - atom name
        "label_alt_id" => (String1, :_print_dot), # nothing (a dot: .)
        "label_comp_id" => (String7, :resname), # redundant - residue name
        "label_asym_id" => (String7, :chain), # redundant - chain id
        "label_entity_id" => (Int32, :_entity_id), # entity_id: print 1 for protein, 2 for water, 3 for all others
        "label_seq_id" => (Int32, :resnum), # redundant - residue number
        "pdbx_PDB_ins_code" => (String1, :_print_question_mark), # nothing - print ?
        "Cartn_x" => (Float32, :x),
        "Cartn_y" => (Float32, :y),
        "Cartn_z" => (Float32, :z),
        "occupancy" => (Float32, :occup),
        "B_iso_or_equiv" => (Float32, :beta),
        "pdbx_formal_charge" => (Float32, :charge),
        "auth_seq_id" => (Int32, :resnum), # Standard mmCIF
        "auth_comp_id" => (String7, :resname), # Standard mmCIF
        "auth_asym_id" => (String7, :chain), # Standard mmCIF
        "auth_atom_id" => (String7, :name), # Standard mmCIF
        "pdbx_PDB_model_num" => (Int32, :model),
    )
    return replace_custom_fields!(_atom_symbol_for_cif_field, field_assignment)
end

function _get_mmcif_field(atom, field_name)
    field_name == :_print_atom && return "ATOM"
    field_name == :_print_dot && return "."
    field_name == :_print_question_mark && return "?"
    field_name == :charge && return round(Int, charge(atom))
    if field_name == :_entity_id
        isprotein(atom) && return 1
        iswater(atom) && return 2
        return 3
    end
    if field_name in fieldnames(typeof(atom))
        return getfield(atom, field_name)
    end
    throw(ArgumentError(":$field_name not found for printing mmcif field in atom"))
end

"""
    write_mmcif(filename, atoms::AbstractVector{<:Atom}, [selection]; field_assignment=nothing)

Write a mmCIF file with the atoms in `atoms` to `filename`. The optional `selection` argument is a string or function
that can be used to select a subset of the atoms in `atoms`. For example, `write_mmcif(atoms, "test.cif", "name CA")`.

The optional `field_assignment` argument is a dictionary that can be used to assign custom fields to the mmCIF file.

"""
function write_mmcif(
    filename::AbstractString,
    atoms::AbstractVector{<:Atom},
    selection::String;
    field_assignment::Union{Nothing,Dict{String,Symbol}}=nothing,
)
    query = parse_query(selection)
    write_mmcif(filename, atoms, atom -> apply_query(query, atom); field_assignment)
end

function write_mmcif(
    filename::AbstractString, 
    atoms::AbstractVector{<:Atom},
    selection_function::Function=all;
    field_assignment::Union{Nothing,Dict{String,Symbol}}=nothing,
)
    _cif_fields = _supported_write_cif_fields(field_assignment)
    open(expanduser(filename), "w") do file
        # Header
        println(file, "_software.name PDBTools.jl $(VERSION)")
        println(file, "_loop")
        println(file, "_atom_site.group_PDB")
        # Atom fields
        for (key, _) in _cif_fields
            println(file, "_atom_site.$key")
        end
        index = 0
        for atom in atoms
            !selection_function(atom) && continue
            index += 1
            buff = IOBuffer(; append=true)
            write(buff, "ATOM $(@sprintf("%9i", index))")
            for (_, field) in _cif_fields
                field_type = first(field)
                field_name = last(field)
                if !(field_name == :index_pdb) && (field_type <: Integer)
                    write(buff, "$(@sprintf("%7i", _get_mmcif_field(atom, field_name)))")
                elseif field_type <: AbstractFloat
                    write(buff, "$(@sprintf("%12.5f",_get_mmcif_field(atom, field_name)))")
                elseif field_type == String1
                    write(buff, "$(@sprintf(" %1s ", _get_mmcif_field(atom, field_name)))")
                elseif field_type == String3
                    write(buff, "$(@sprintf(" %3s ", _get_mmcif_field(atom, field_name)))")
                elseif field_type == String7
                    write(buff, "$(@sprintf(" %7s ", _get_mmcif_field(atom, field_name)))")
                end
            end
            println(file, String(take!(buff)))
        end
        println(file, "#")
    end
end

@testitem "write_mmcif" begin
    using PDBTools
    ats = read_pdb(PDBTools.SMALLPDB)
    tmpfile = tempname() * ".cif"
    write_mmcif(tmpfile, ats)
    @test isfile(tmpfile)
    ats_cif = read_mmcif(tmpfile)
    @test all(position(at1) ≈ position(at2) for (at1, at2) in zip(ats, ats_cif))
    field_assignment = Dict("test" => :name)
    write_mmcif(tmpfile, ats; field_assignment)
    ats0 = read_mmcif(tmpfile)
    @test all(name(at) == "X" for at in ats0)
    ats1 = read_mmcif(tmpfile; field_assignment)
    @test all(name(at1) == name(at2) for (at1, at2) in zip(ats, ats1))
    field_assignment = Dict("test" => :abc)
    @test_throws ArgumentError write_mmcif(tmpfile, ats; field_assignment)

    cA = select(ats, at -> name(at) == "CA")
    write_mmcif(tmpfile, ats, "name CA")
    @test isfile(tmpfile)
    ats_cif = read_mmcif(tmpfile)
    @test all(position(at1) ≈ position(at2) for (at1, at2) in zip(cA, ats_cif))

    write_mmcif(tmpfile, ats, at -> name(at) == "CA")
    @test isfile(tmpfile)
    ats_cif = read_mmcif(tmpfile)
    @test all(position(at1) ≈ position(at2) for (at1, at2) in zip(cA, ats_cif))
end
