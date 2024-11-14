"""
    write_mmcif(filename, atoms::AbstractVector{<:Atom}, [selection])

Write a mmCIF file with the atoms in `atoms` to `filename`. The optional `selection` argument is a string
that can be used to select a subset of the atoms in `atoms`. For example, `write_mmcif(atoms, "test.cif", "name CA")`.

"""
function write_mmcif(filename::AbstractString, atoms::AbstractVector{<:Atom}, selection::String)
    query = parse_query(selection)
    write_mmcif(filename, atoms; only=atom -> apply_query(query, atom))
end

function write_mmcif(filename::AbstractString, atoms::AbstractVector{<:Atom}; only::Function=all)
    _cif_fields = _supported_cif_fields()
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
            !only(atom) && continue
            index += 1
            buff = IOBuffer(;append=true)
            write(buff, "ATOM $(@sprintf("%9i", index))")
            for (_, field) in _cif_fields
                field_type = first(field)
                field_name = last(field)
                if !(field_name == :index_pdb) && (field_type <: Integer)
                    write(buff, "$(@sprintf("%7i", getfield(atom, field_name)))")
                elseif field_type <: AbstractFloat
                    write(buff, "$(@sprintf("%12.5f", getfield(atom, field_name)))")
                elseif field_type <: AbstractString
                    write(buff, "$(@sprintf("%7s", getfield(atom, field_name)))")
                end
            end
            println(file, String(take!(buff)))
        end
    end
end

@testitem "write_mmcif" begin
    using PDBTools
    ats = read_pdb(PDBTools.SMALLPDB)
    tmpfile = tempname()*".cif"
    write_mmcif(tmpfile, ats)
    @test isfile(tmpfile)
    ats_cif = read_mmcif(tmpfile)
    @test all(position(at1) ≈ position(at2) for (at1, at2) in zip(ats, ats_cif))
end
