"""
    read_mmcif(mmCIF_file::String, selection::String)
    read_mmcif(mmCIF_file::String; only::Function = all)

    read_mmcif(mmCIF_data::IOBuffer, selection::String)
    read_mmcif(mmCIF_data::IOBuffer; only::Function = all)

Reads a mmCIF file and stores the data in a vector of type `Atom`. 

If a selection is provided, only the atoms matching the selection will be read. 
For example, `resname ALA` will select all the atoms in the residue ALA.

If the `only` function keyword is provided, only the atoms for which `only(atom)` is true will be returned.

### Examples

```julia-repl
julia> ats = read_mmcif(PDBTools.SMALLCIF)
   Array{Atoms,1} with 7 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     VAL     A        1        1    6.204   16.869    4.854  1.00 49.05     1                 1
       2   CA     VAL     A        1        1    6.913   17.759    4.607  1.00 43.14     1                 2
       3    C     VAL     A        1        1    8.504   17.378    4.797  1.00 24.80     1                 3
       5   CB     VAL     A        1        1    6.369   19.044    5.810  1.00 72.12     1                 5
       6  CG1     VAL     A        1        1    7.009   20.127    5.418  1.00 61.79     1                 6
       7  CG2     VAL     A        1        1    5.246   18.533    5.681  1.00 80.12     1                 7

julia> ats = read_mmcif(PDBTools.SMALLCIF, "index < 3")
   Array{Atoms,1} with 2 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     VAL     A        1        1    6.204   16.869    4.854  1.00 49.05     1                 1
       2   CA     VAL     A        1        1    6.913   17.759    4.607  1.00 43.14     1                 2

julia> ats = read_mmcif(PDBTools.SMALLCIF; only = at -> name(at) == "CA")
   Array{Atoms,1} with 1 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       2   CA     VAL     A        1        1    6.913   17.759    4.607  1.00 43.14     1                 2

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
    if mod(length(atoms), 1000) == 0
        if Sys.free_memory() < (1 - memory_available) * Sys.total_memory()
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
function _supported_cif_fields()
    # We need this to be indexable (not a Dict) to keep the order of the fields
    # when writting the mmCIF file
    return (
        "id" => (Int32, :index_pdb),
        "label_atom_id" => (String7, :name), # VMD
        "label_comp_id" => (String7, :resname), # VMD
        "label_asym_id" => (String3, :chain), # VMD
        "label_seq_id" => (Int32, :resnum), # VMD
        "Cartn_x" => (Float32,:x),
        "Cartn_y" => (Float32,:y),
        "Cartn_z" => (Float32,:z),
        "occupancy" => (Float32,:occup),
        "B_iso_or_equiv" => (Float32,:beta),
        "pdbx_formal_charge" => (Float32,:charge),
        "pdbx_PDB_model_num" => (Int32,:model),
        "auth_atom_id" => (String7, :name), # Standard mmCIF
        "auth_comp_id" => (String7, :resname), # Standard mmCIF
        "auth_asym_id" => (String3, :chain), # Standard mmCIF
        "auth_seq_id" => (Int32, :resnum), # Standard mmCIF
    )
end

function _parse_mmCIF(
    cifdata::Union{IOStream,IOBuffer};
    only::Function,
    memory_available::Real=1, # defaults to 100% (MacOS requires this)
    stop_at=nothing,
)
    _atom_symbol_for_cif_field = _supported_cif_fields()
    _atom_site_field_inds = Dict{String,Int}()
    ifield = 0
    atoms = Atom{Nothing}[]
    lastatom = Atom{Nothing}()
    _atom_field_columns = Vector{Tuple{Int,Tuple{DataType,Symbol}}}()
    local NCOLS, col_indices, col_field
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
    b = @benchmark read_mmcif($(PDBTools.SMALLCIF)) samples=1 evals=1
    @test b.allocs < 200
    ats = read_mmcif(PDBTools.TESTCIF)
    @test count(at -> resname(at) == "HOH", ats) == 445
    @test count(at -> isprotein(at), ats) == 2966
    @test length(eachresidue(select(ats, "protein"))) == 354
    tmpfile = tempname()*".cif"
    write_mmcif(tmpfile, ats, "protein")
    prot = read_mmcif(tmpfile)
    @test length(prot) == 2966
    ats = read_mmcif(PDBTools.TESTCIF, "resname HOH")
    @test length(ats) == 445

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
end