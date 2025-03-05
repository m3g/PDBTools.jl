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
       3    C     SER     A      225        1   45.257   81.872   70.236  1.00 67.90     1                 3
       4    O     SER     A      225        1   45.823   80.796   69.974  1.00 64.85     1                 4
       5   CB     SER     A      225        1   47.147   82.980   71.413  1.00 70.79     1                 5
       6   OG     SER     A      225        1   46.541   82.639   72.662  1.00 73.55     1                 6
⋮
    1863    N     THR     A      462      238  -26.699   75.053   51.551  1.00 81.01     1              1863
    1864   CA     THR     A      462      238  -26.915   73.645   51.198  1.00 82.96     1              1864
    1865    C     THR     A      462      238  -26.026   72.767   52.099  1.00 85.11     1              1865
    1866    O     THR     A      462      238  -25.993   73.062   53.309  1.00 88.20     1              1866
    1867   CB     THR     A      462      238  -26.675   73.396   49.649  1.00 81.31     1              1867
    1868  OG1     THR     A      462      238  -27.462   74.325   48.885  1.00 79.98     1              1868
    1869  CG2     THR     A      462      238  -27.063   71.965   49.222  1.00 78.62     1              1869
    1870  OXT     THR     A      462      238  -25.379   71.816   51.613  1.00 84.35     1              1870

```
"""
function wget(pdb_id::String, selection::String; format::AbstractString="mmCIF")
    query = parse_query(selection)
    return wget(pdb_id, only=atom -> apply_query(query, atom); format)
end

function _wget(pdb_id, format, only)
    buf = IOBuffer()
    atoms = try
        Downloads.download("https://files.rcsb.org/download/$(uppercase(pdb_id)).$format", buf)
        seekstart(buf)
        if format == "pdb"
            read_pdb(buf, only=only)
        else
            read_mmcif(buf, only=only)
        end
    catch
        @info "Failed downloading from `download` PDB repository, trying `view` repository ..."
        Downloads.download("https://files.rcsb.org/view/$(uppercase(pdb_id)).$format", buf)
        seekstart(buf)
        if format == "pdb"
            read_pdb(buf, only=only)
        else
            read_mmcif(buf, only=only)
        end
    end
    return atoms
end


function wget(pdb_id::String; only=all, format::Union{AbstractString,Nothing}=nothing)
    atoms = if format == "PDB"
        _wget(pdb_id, "pdb", only)
    elseif isnothing(format) || format == "mmCIF" 
        _wget(pdb_id, "cif", only)
    else
        throw(ArgumentError("""\n
            format must be either "PDB" or "mmCIF"
        
        """))
    end
    return atoms
end

@testitem "wget" begin
    using PDBTools
    protein = wget("1LBD","chain A")
    @test length(protein) == 1870
    protein = wget("1LBD","chain A"; format="PDB")
    @test length(protein) == 1870
    protein = wget("1LBD","chain A"; format="mmCIF")
    @test length(protein) == 1870
    @test_throws ArgumentError wget("1LBD","chain A"; format="mmcif")
end
