"""
    wget(PDBid; selection; format="mmCIF")

Retrieves a PDB file from the protein data bank. Selections may be applied.

The optional format argument can be either "mmCIF" or "PDB". The default is "mmCIF".
To download the data of large structures, it is recommended to use the "mmCIF" format.

### Example

```jldoctest
julia> using PDBTools

julia> protein = wget("1LBD","chain A")
   Vector{Atom{Nothing}} with 1870 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     SER     A      225        1   45.228   84.358   70.638  1.00 67.05     1                 1
       2   CA     SER     A      225        1   46.080   83.165   70.327  1.00 68.73     1                 2
⋮
    1869  CG2     THR     A      462      238  -27.063   71.965   49.222  1.00 78.62     1              1869
    1870  OXT     THR     A      462      238  -25.379   71.816   51.613  1.00 84.35     1              1870

```
"""
function wget(pdb_id::String, selection::String; format::AbstractString="mmCIF")
    query = parse_query(selection)
    return wget(pdb_id, atom -> apply_query(query, atom); format)
end

function _wget(pdb_id, selection_function::Function; format)
    buf = IOBuffer()
    atoms = try
        Downloads.download("https://files.rcsb.org/download/$(uppercase(pdb_id)).$format", buf)
        seekstart(buf)
        if format == "pdb"
            read_pdb(buf, selection_function)
        else
            read_mmcif(buf, selection_function)
        end
    catch
        @info "Failed downloading from `download` PDB repository, trying `view` repository ..."
        Downloads.download("https://files.rcsb.org/view/$(uppercase(pdb_id)).$format", buf)
        seekstart(buf)
        if format == "pdb"
            read_pdb(buf, selection_function)
        else
            read_mmcif(buf, selection_function)
        end
    end
    return atoms
end


function wget(pdb_id::String, selection_function::Function=all; format::Union{AbstractString,Nothing}=nothing)
    atoms = if format == "PDB"
        _wget(pdb_id, selection_function; format="pdb")
    elseif isnothing(format) || format == "mmCIF"
        _wget(pdb_id, selection_function; format="cif")
    else
        throw(ArgumentError("""\n
            format must be either "PDB" or "mmCIF"

        """))
    end
    return atoms
end

@testitem "wget" begin
    using PDBTools
    protein = wget("1LBD", "chain A")
    @test length(protein) == 1870
    protein = wget("1LBD", "chain A"; format="PDB")
    @test length(protein) == 1870
    protein = wget("1LBD", "chain A"; format="mmCIF")
    @test length(protein) == 1870
    protein = wget("1LBD", at -> chain(at) == "A"; format="mmCIF")
    @test length(protein) == 1870
    @test_throws ArgumentError wget("1LBD", "chain A"; format="mmcif")
end
