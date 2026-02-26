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
        "type_symbol" => (StringType, :pdb_element), # Standard mmCIF
        "label_atom_id" => (StringType, :name), # redundant - atom name
        "label_alt_id" => (String1, :_print_dot), # nothing (a dot: .)
        "label_comp_id" => (StringType, :resname), # redundant - residue name
        "label_asym_id" => (StringType, :chain), # redundant - chain id
        "label_entity_id" => (Int32, :_entity_id), # entity_id: print 1 for protein, 2 for water, 3 for all others
        "label_seq_id" => (Int32, :resnum), # redundant - residue number
        "pdbx_PDB_ins_code" => (String1, :_print_question_mark), # nothing - print ?
        "Cartn_x" => (Float32, :x),
        "Cartn_y" => (Float32, :y),
        "Cartn_z" => (Float32, :z),
        "occupancy" => (Float32, :occup),
        "B_iso_or_equiv" => (Float32, :beta),
        "pdbx_formal_charge" => (String1, :_print_question_mark), # use ? for unknown/zero charges
        "auth_seq_id" => (Int32, :resnum), # Standard mmCIF
        "auth_comp_id" => (StringType, :resname), # Standard mmCIF
        "auth_asym_id" => (StringType, :chain), # Standard mmCIF
        "auth_atom_id" => (StringType, :name), # Standard mmCIF
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

const _fmt_string = generate_formatter("%$(_length(StringType))s ")
const _fmt_int = generate_formatter("%9i")
const _fmt_float = generate_formatter("%12.5f")
_fmt(::Type{<:InlineString}, val) = _fmt_string(val)
_fmt(::Type{<:Integer}, val) = _fmt_int(val)
_fmt(::Type{<:AbstractFloat}, val) = _fmt_float(val)

function write_mmcif(
    filename::AbstractString, 
    atoms::AbstractVector{<:Atom},
    selection_function::Function=all;
    field_assignment::Union{Nothing,Dict{String,Symbol}}=nothing,
)
    _cif_fields = _supported_write_cif_fields(field_assignment)
    open(expanduser(filename), "w") do file
        # Data block header
        println(file, "data_MODEL")
        println(file, "#")
        println(file, "_entry.id MODEL")
        println(file, "#")
        println(file, "_audit_conform.dict_name       mmcif_pdbx.dic")
        println(file, "_audit_conform.dict_version    5.397")
        println(file, "_audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic")
        println(file, "#")
        # Collect chain→entity_id (unique per chain) and unique (chain, resnum) pairs
        chain_entity = OrderedDict{String,Int}()
        chain_type = OrderedDict{String,Symbol}()  # :protein, :water, :other
        residue_list = OrderedDict{Tuple{String,Int32},Nothing}()
        for atom in atoms
            !selection_function(atom) && continue
            ch = string(atom.chain)
            if !haskey(chain_entity, ch)
                chain_entity[ch] = length(chain_entity) + 1
                chain_type[ch] = isprotein(atom) ? :protein : iswater(atom) ? :water : :other
            end
            key = (ch, atom.resnum)
            haskey(residue_list, key) || (residue_list[key] = nothing)
        end
        # _entity: one entry per chain (each chain gets its own entity_id to avoid DSSP hang)
        println(file, "loop_")
        println(file, "_entity.id")
        println(file, "_entity.type")
        for (ch, eid) in chain_entity
            type_str = chain_type[ch] == :protein ? "polymer" : chain_type[ch] == :water ? "water" : "non-polymer"
            println(file, "$eid $type_str")
        end
        println(file, "#")
        # _entity_poly (only for polymer/protein chains)
        if any(t -> t == :protein, values(chain_type))
            println(file, "loop_")
            println(file, "_entity_poly.entity_id")
            println(file, "_entity_poly.type")
            for (ch, eid) in chain_entity
                chain_type[ch] == :protein && println(file, "$eid polypeptide(L)")
            end
            println(file, "#")
        end
        # _struct_asym
        println(file, "loop_")
        println(file, "_struct_asym.id")
        println(file, "_struct_asym.entity_id")
        for (ch, eid) in chain_entity
            println(file, "$ch $eid")
        end
        println(file, "#")
        # _pdbx_poly_seq_scheme (only for protein/polymer chains)
        poly_residues = [(ch, rnum) for (ch, rnum) in keys(residue_list) if get(chain_type, ch, :other) == :protein]
        if !isempty(poly_residues)
            println(file, "loop_")
            println(file, "_pdbx_poly_seq_scheme.asym_id")
            println(file, "_pdbx_poly_seq_scheme.seq_id")
            println(file, "_pdbx_poly_seq_scheme.pdb_strand_id")
            println(file, "_pdbx_poly_seq_scheme.pdb_seq_num")
            println(file, "_pdbx_poly_seq_scheme.pdb_ins_code")
            for (ch, rnum) in poly_residues
                println(file, "$ch $rnum $ch $rnum ?")
            end
            println(file, "#")
        end
        # _atom_site loop
        println(file, "loop_")
        println(file, "_atom_site.group_PDB")
        for (key, _) in _cif_fields
            println(file, "_atom_site.$key")
        end
        index = 0
        for atom in atoms
            !selection_function(atom) && continue
            index += 1
            buff = IOBuffer(; append=true)
            write(buff, "ATOM $(_fmt(Int, index))")
            for (_, field) in _cif_fields
                field_type = first(field)
                field_name = last(field)
                field_name == :index_pdb && continue
                val = field_name == :_entity_id ? chain_entity[string(atom.chain)] : _get_mmcif_field(atom, field_name)
                write(buff, _fmt(field_type, val))
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
