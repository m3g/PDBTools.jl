"""
    read_mmcif(mmCIF_file::String, selection::String; field_assignment)
    read_mmcif(mmCIF_file::String; only::Function = all, field_assignment)

    read_mmcif(mmCIF_data::IOBuffer, selection::String; field_assignment)
    read_mmcif(mmCIF_data::IOBuffer; only::Function = all, field_assignment)

Reads a mmCIF file and stores the data in a vector of type `Atom`. 

All fields except the file name are optional.

If a selection is provided, only the atoms matching the selection will be read. 
For example, `resname ALA` will select all the atoms in the residue ALA.

If the `only` function keyword is provided, only the atoms for which `only(atom)` is true will be returned.

The `field_assignment` keyword is `nothing` (default) or a `Dict{String,Symbol}` and can be used to specify which fields in the mmCIF file should be read into the `Atom` type.
For example `field_assignment = Dict("type_symbol" => :name)` will read the `_atom_site.type_symbol` field in the mmCIF 
file into the `name` field of the `Atom` type.

The default assignment is follows the standard mmCIF convention:

```julia
Dict{String,Symbol}(
    "id" => :index_pdb
    "Cartn_x" => :x
    "Cartn_y" => :y
    "Cartn_z" => :z
    "occupancy" => :occup
    "B_iso_or_equiv" => :beta
    "pdbx_formal_charge" => :charge
    "pdbx_PDB_model_num" => :model
    "label_atom_id" => :name
    "label_comp_id" => :resname
    "label_asym_id" => :chain
    "auth_seq_id" => :resnum
    "type_symbol" => :pdb_element
)
```

Source: https://mmcif.wwpdb.org/docs/tutorials/content/atomic-description.html

### Examples

```jldoctest
julia> using PDBTools

julia> ats = read_mmcif(PDBTools.TESTCIF)
   Vector{Atom{Nothing}} with 76 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     GLY     A        1        1   -4.564   25.503   24.113  1.00 24.33     1                 1
       2   CA     GLY     A        1        1   -4.990   26.813   24.706  1.00 24.29     1                 2
       3    C     GLY     A        1        1   -4.558   27.997   23.861  1.00 23.83     1                 3
                                                       ⋮
      74    O     HOH     Q       62       14    2.156   31.115   23.421  1.00 18.43     1              2979
      75    O     HOH     Q       63       15   -3.585   34.725   20.903  1.00 19.82     1              2980
      76    O     HOH     Q       64       16   -4.799   40.689   37.419  1.00 20.13     1              2981

julia> ats = read_mmcif(PDBTools.TESTCIF, "index < 3")
   Vector{Atom{Nothing}} with 2 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     GLY     A        1        1   -4.564   25.503   24.113  1.00 24.33     1                 1
       2   CA     GLY     A        1        1   -4.990   26.813   24.706  1.00 24.29     1                 2

julia> ats = read_mmcif(PDBTools.TESTCIF; only = at -> name(at) == "CA")
   Vector{Atom{Nothing}} with 11 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       2   CA     GLY     A        1        1   -4.990   26.813   24.706  1.00 24.29     1                 2
       6   CA     GLN     A        2        2   -4.738   30.402   23.484  1.00 23.74     1                 6
      15   CA     ARG     A        3        3   -1.234   31.761   22.771  1.00 23.94     1                15
                                                       ⋮
      63   CA     ASP     A        9        9    7.831   38.316   36.914  1.00 24.27     1                63
      70   CA      CA     G     1003       10  -24.170   27.201   64.364  1.00 27.40     1              2967
      71   CA      CA     H     1004       11  -10.624   32.854   69.292  1.00 29.53     1              2968

```

"""
function read_mmcif end

function read_mmcif(file::Union{String,IOBuffer}, selection::String; kargs...)
    query = parse_query(selection)
    return read_mmcif(file, only=atom -> apply_query(query, atom); kargs...)
end

function read_mmcif(cifdata::IOBuffer; only::Function=all, kargs...)
    atoms = _parse_mmCIF(cifdata; only, kargs...)
    return atoms
end

function read_mmcif(file::String; only=all, kargs...)
    atoms = open(expanduser(file), "r") do f
        _parse_mmCIF(f; only, kargs...)
    end
    return atoms
end

function _maximum_read(atoms, stop_at, memory_available)
    if !isnothing(stop_at) && (length(atoms) >= stop_at) 
        return true
    end
    if length(atoms) > 0 && mod(length(atoms), 50) == 0
        estimated_bytes = length(atoms) * Base.summarysize(atoms[1])
        if estimated_bytes > memory_available * Sys.total_memory()
            @warn """\n
                Memory limit reached. $(length(atoms)) atoms read so far will be returned.
                Size of the atoms array: $(round(Base.summarysize(atoms) / 1024^2; digits=3)) MB

            """ _file = nothing _line = nothing
            return true
        end
    end
    return false
end

function _error_if_no_atoms(atoms)
    if length(atoms) == 0
        throw(ArgumentError("""\n 
            Could not find any atom in mmCIF file matching the selection. 

        """))
    end
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
function _supported_cif_fields(field_assignment)
    # We need this to be indexable (not a Dict) to keep the order of the fields
    # when writing the mmCIF file
    _atom_symbol_for_cif_field = OrderedDict{String, Tuple{DataType, Symbol}}(
        "id" => (Int32, :index_pdb),
        "Cartn_x" => (Float32,:x),
        "Cartn_y" => (Float32,:y),
        "Cartn_z" => (Float32,:z),
        "occupancy" => (Float32,:occup),
        "B_iso_or_equiv" => (Float32,:beta),
        "pdbx_formal_charge" => (Float32,:charge),
        "pdbx_PDB_model_num" => (Int32,:model),
        "label_atom_id" => (String7, :name), # Standard mmCIF
        "label_comp_id" => (String7, :resname), # Standard mmCIF
        "label_asym_id" => (String3, :chain), # Standard mmCIF
        "auth_seq_id" => (Int32, :resnum), # Standard mmCIF
        "type_symbol" => (String7, :pdb_element), # Standard mmCIF
        #"auth_atom_id" => (String7, :name), # alternate - standard mmCIF
        #"auth_comp_id" => (String7, :resname), # alternate - standard mmCIF
        #"auth_asym_id" => (String3, :chain), # alternate = standard mmCIF
    )
    if !isnothing(field_assignment)
        for (custom_key, custom_field) in field_assignment 
            if !(custom_field in fieldnames(Atom))
                throw(ArgumentError("""\n
                    Field $custom_field not available in the PDBTools.Atom type.

                """))
            end
            if custom_field == :custom
                throw(ArgumentError("""\n
                    Setting a custom field to :custom is not currently supported.

                """))
            end
            # Type to assign to the field
            T = fieldtype(Atom, custom_field)
            # Remove the entry from the standard mmCIF fields, if available
            if haskey(_atom_symbol_for_cif_field, custom_key)
                pop!(_atom_symbol_for_cif_field, custom_key)
            end
            # Remove entry corresponding to custom_field, if available
            for p in _atom_symbol_for_cif_field
                key, val = p
                if last(val) == custom_field
                    pop!(_atom_symbol_for_cif_field, key)
                end
            end
            # Add entry associated to new assignment
            push!(_atom_symbol_for_cif_field, custom_key => (T, custom_field))
        end
    end
    return _atom_symbol_for_cif_field
end

@testitem "_supported_cif_fields" begin
    using PDBTools 
    using InlineStrings
    field_assignment = Dict(
        "type_symbol" => :name, 
        "label_comp_id" => :resname,
        "label_asym_id" => :chain, 
        "label_seq_id" => :resnum, 
    )
    _cif_fields = PDBTools._supported_cif_fields(field_assignment)
    @test length(_cif_fields) == 12
    @test _cif_fields["id"] == (Int32, :index_pdb)
    @test _cif_fields["type_symbol"] == (String7, :name)
    @test _cif_fields["label_comp_id"] == (String7, :resname)
    @test _cif_fields["label_asym_id"] == (String3, :chain)
    @test _cif_fields["label_seq_id"] == (Int32, :resnum)
    @test _cif_fields["Cartn_x"] == (Float32, :x)
    @test _cif_fields["Cartn_y"] == (Float32, :y)
    @test _cif_fields["Cartn_z"] == (Float32, :z)
    @test _cif_fields["occupancy"] == (Float32, :occup)
    @test _cif_fields["B_iso_or_equiv"] == (Float32, :beta)
    @test _cif_fields["pdbx_formal_charge"] == (Float32, :charge)
    @test _cif_fields["pdbx_PDB_model_num"] == (Int32, :model)
    push!(field_assignment, "type_symbol" => :fail)
    @test_throws ArgumentError PDBTools._supported_cif_fields(field_assignment)
    field_assignment = Dict("type_symbol" => :custom)
    @test_throws ArgumentError PDBTools._supported_cif_fields(field_assignment)
end

function _parse_mmCIF(
    cifdata::Union{IOStream,IOBuffer};
    only::Function,
    memory_available::Real=0.8,
    stop_at=nothing,
    field_assignment::Union{Nothing,Dict{String,Symbol}} = nothing,
)
    _atom_symbol_for_cif_field = _supported_cif_fields(field_assignment)
    atoms = Atom{Nothing}[]
    lastatom = Atom{Nothing}()
    _atom_field_columns = Vector{Tuple{Int,Tuple{DataType,Symbol}}}()
    local NCOLS, col_indices, col_field
    ifield = 0
    _atom_site_field_inds = Dict{String,Int}()
    for line in eachline(cifdata)
        # Reading the headers of the _atom_site loop
        if occursin("loop_", line)
            ifield = 0
            empty!(_atom_site_field_inds)
        end
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
            for (key, keyval) in _atom_symbol_for_cif_field
                if haskey(_atom_site_field_inds, key)
                    push!(_atom_field_columns, (_atom_site_field_inds[key], keyval))
                end
            end
            sort!( _atom_field_columns; by = first)
            col_indices = NTuple{length(_atom_field_columns),Int}(first(el) for el in _atom_field_columns)
            col_field = NTuple{length(_atom_field_columns),Symbol}(last(last(el)) for el in _atom_field_columns)
            NCOLS = length(keys(_atom_site_field_inds))
            break
        end
    end
    inds_and_names = ntuple(length(col_indices)) do i 
        ((col_indices[i], Val(col_field[i]))) 
    end
    seekstart(cifdata)
    for line in eachline(cifdata)
        if startswith(line, r"ATOM|HETATM")
            atom = read_atom_mmcif(Val(NCOLS), line, inds_and_names, lastatom)
            only(atom) && push!(atoms, atom)
            _maximum_read(atoms, stop_at, memory_available) && break
            lastatom = atom
        end
    end
    seekstart(cifdata)
    _error_if_no_atoms(atoms)
    return atoms
end

function read_atom_mmcif(::Val{NCOLS}, record, inds_and_names, lastatom::Atom) where {NCOLS}
    field_values = NTuple{NCOLS}(eachsplit(record))
    atom = Atom{Nothing}(; index = index(lastatom) + 1, residue = residue(lastatom))
    _fast_setfield!(atom, field_values, inds_and_names)
    if !same_residue(atom, lastatom)
        atom.residue = residue(lastatom) + 1
    end
    return atom
end

@testitem "read_mmcif" begin
    using PDBTools
    using BenchmarkTools
    b = @benchmark read_mmcif($(PDBTools.TESTCIF)) samples=1 evals=1
    @test b.allocs < 500
    ats = read_mmcif(PDBTools.TESTCIF)
    @test count(iswater, ats) == 5
    @test count(isprotein, ats) == 69
    @test length(eachresidue(filter(isprotein, ats))) == 9
    tmpfile = tempname()*".cif"
    write_mmcif(tmpfile, ats, "protein")
    prot = read_mmcif(tmpfile)
    @test length(prot) == 69
    ats = read_mmcif(PDBTools.TESTCIF, "resname HOH")
    @test length(ats) == 5

    # peformance tests for innner functions for reading cif fields into Atoms
    record = "ATOM   1    N  N   . VAL A 1 1   ? 6.204   16.869  4.854   1.00 49.05 ? 1   VAL A N   1"
    inds_and_names = ((2, Val{:index_pdb}()), (4, Val{:name}()), (6, Val{:resname}()), (7, Val{:chain}()), (9, Val{:resnum}()), (11, Val{:x}()), (12, Val{:y}()), (13, Val{:z}()), (14, Val{:occup}()), (15, Val{:beta}()), (16, Val{:charge}()), (17, Val{:resnum}()), (18, Val{:resname}()), (19, Val{:chain}()), (20, Val{:name}()), (21, Val{:model}()))
    lastatom = Atom()
    NCOLS = 21
    b = @benchmark PDBTools.read_atom_mmcif($(Val(NCOLS)), $record, $inds_and_names, $lastatom) samples=1 evals=1
    @test b.allocs == 1
    field_values = NTuple{NCOLS}(eachsplit(record))
    atom = Atom{Nothing}(; index = index(lastatom) + 1, residue = residue(lastatom))
    b = @benchmark PDBTools._fast_setfield!($atom, $field_values, $inds_and_names) samples=1 evals=1
    @test b.allocs == 0

    # Test early stoppers
    ats = read_mmcif(PDBTools.TESTCIF; stop_at=10)
    @test length(ats) == 10
    ats = read_mmcif(PDBTools.TESTCIF; memory_available=1e-10)
    @test length(ats) == 50 

    # Read custom fields instead
    field_assignment = Dict(
        "type_symbol" => :name, 
        "label_comp_id" => :resname, 
        "label_asym_id" => :chain, 
        "label_seq_id" => :resnum, 
    )
    ats0 = read_mmcif(PDBTools.TESTCIF)
    ats1 = read_mmcif(PDBTools.TESTCIF; field_assignment)
    @test all(name.(ats1) .== pdb_element.(ats0))
    @test all(resname.(ats1) .== resname.(ats0))
    @test all(chain.(ats1) .== chain.(ats0))
    @test all(resnum.(filter(isprotein,ats1)) .== resnum.(filter(isprotein,ats0)))
    @test resnum.(filter(iswater, ats0)) == Int32[60, 61, 62, 63, 64]
    @test resnum.(filter(iswater, ats1)) == Int32[0, 0, 0, 0, 0] 
end