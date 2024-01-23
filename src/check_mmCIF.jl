#
# Structure to read mmCIF field data indices
#
Base.@kwdef mutable struct indices_mmCIF_fields
    index::Int = 0
    name::Int = 0
    resname::Int = 0
    chain::Int = 0
    resnum::Int = 0
    x::Int = 0
    y::Int = 0
    z::Int = 0
    beta::Int = 0
    occup::Int = 0
end

#
# Check if structure is in mmCIF format and returns the list of 
# desired fields if that is the case
#
function check_mmCIF(file::String) 
    data = open(file, "r")
    mmCIF, mmCIF_fields = _check_mmCIF(data)
    close(data)
    return mmCIF, mmCIF_fields
end
check_mmCIF(data::IOBuffer) = _check_mmCIF(data)

function _check_mmCIF(data::Union{IOStream, IOBuffer})
    local ifield, mmCIF, line
    mmCIF_fields = indices_mmCIF_fields()
    mmCIF = false
    for line_string in eachline(data)
        line = strip(line_string)
        if length(line) < 5
            continue
        end
        if line[1:5] == "loop_"
            mmCIF = true
            ifield = 0
            while !eof(data)
                try
                    line = strip(readline(data))
                catch
                    break
                end
                if line[1:1] == "#"
                    continue
                end
                if ifield > 0
                    ifield = ifield + 1
                    if occursin("_atom_site.label_atom_id", line)
                        mmCIF_fields.name = ifield
                    end
                    if occursin("_atom_site.id", line)
                        mmCIF_fields.index = ifield
                    end
                    if occursin("_atom_site.label_comp_id", line)
                        mmCIF_fields.resname = ifield
                    end
                    if occursin("_atom_site.label_asym_id", line)
                        mmCIF_fields.chain = ifield
                    end
                    if occursin("_atom_site.label_seq_id", line)
                        mmCIF_fields.resnum = ifield
                    end
                    if occursin("_atom_site.Cartn_x", line)
                        mmCIF_fields.x = ifield
                    end
                    if occursin("_atom_site.Cartn_y", line)
                        mmCIF_fields.y = ifield
                    end
                    if occursin("_atom_site.Cartn_z", line)
                        mmCIF_fields.z = ifield
                    end
                    if occursin("_atom_site.B_iso_or_equiv", line)
                        mmCIF_fields.beta = ifield
                    end
                    if occursin("_atom_site.occupancy", line)
                        mmCIF_fields.occup = ifield
                    end
                    if occursin("ATOM", line) && occursin("HETATOM", line)
                        break
                    end
                end
                if occursin("_atom_site.group_PDB", line)
                    ifield = 1
                end
            end
        end
    end
    seekstart(data)
    if mmCIF
        if mmCIF_fields.name == 0 ||
           mmCIF_fields.resname == 0 ||
           mmCIF_fields.chain == 0 ||
           mmCIF_fields.resnum == 0 ||
           mmCIF_fields.x == 0 ||
           mmCIF_fields.y == 0 ||
           mmCIF_fields.z == 0
            error(" ERROR: Could not find all necessary fields in mmCIF file. ")
        end
    end
    return mmCIF, mmCIF_fields
end

@testitem "check_mmCIF" begin
    pdbfile = "$(@__DIR__)/../test/structure.pdb"
    @test PDBTools.check_mmCIF(pdbfile)[1] == false
    pdb_data = read(pdbfile, String)
    @test PDBTools.check_mmCIF(IOBuffer(pdb_data))[1] == false
end
