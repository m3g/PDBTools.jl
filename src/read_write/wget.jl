"""
    wget(PDBid; selection; format="mmCIF", assembly=nothing)

Retrieves a PDB file from the protein data bank. Selections may be applied.

The optional format argument can be either "mmCIF" or "PDB". The default is "mmCIF".
To download the data of large structures, it is recommended to use the "mmCIF" format.

By default, the asymmetric unit is downloaded. To retrieve a biological assembly
instead (e.g. the crystallographic tetramer, when the asymmetric unit only contains
one subunit), pass its number as the `assembly` keyword, e.g. `assembly=1`.

RCSB numbers the biological assemblies it has on file for each entry (`assembly=1`,
`assembly=2`, ...); an entry may have none, one, or several, since some depositions
include more than one proposed oligomeric state (for example, a dimer as assembly 1
and a hexamer as assembly 2). `assembly` therefore selects *which* deposited assembly
to fetch, not whether to fetch one, so it takes an assembly number rather than a
`Bool`. Passing a number for an assembly the entry doesn't have raises an error from
the underlying download (RCSB returns a 404).

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
function wget(pdb_id::AbstractString, selection::AbstractString; format::AbstractString="mmCIF", assembly::Union{Integer,Nothing}=nothing)
    return wget(pdb_id, parse_query(selection); format, assembly)
end

# The legacy PDB format names assembly files as "$(id).pdb$(n)" (no dash), while
# mmCIF names them "$(id)-assembly$(n).cif": both confirmed against files.rcsb.org.
function _assembly_url(repository, pdb_id, format, assembly)
    id = uppercase(pdb_id)
    filename = isnothing(assembly) ? "$id.$format" : format == "pdb" ? "$id.pdb$assembly" : "$id-assembly$assembly.$format"
    return "https://files.rcsb.org/$repository/$filename"
end

function _wget(pdb_id, selection_function::Function; format, assembly=nothing)
    buf = IOBuffer()
    atoms = try
        Downloads.download(_assembly_url("download", pdb_id, format, assembly), buf)
        seekstart(buf)
        if format == "pdb"
            read_pdb(buf, selection_function)
        else
            read_mmcif(buf, selection_function)
        end
    catch
        try 
            @info "Failed downloading from `download` PDB repository, trying `view` repository ..."
            Downloads.download(_assembly_url("view", pdb_id, format, assembly), buf)
            seekstart(buf)
            if format == "pdb"
                read_pdb(buf, selection_function)
            else
                read_mmcif(buf, selection_function)
            end
        catch
            error("""\n
                Failed downloading $pdb_id - please check the identifier code or internet connection.
            
            """)
        end
    end
    return atoms
end


function wget(pdb_id::AbstractString, selection_function::Function=all; format::Union{AbstractString,Nothing}=nothing, assembly::Union{Integer,Nothing}=nothing)
    atoms = if format == "PDB"
        _wget(pdb_id, selection_function; format="pdb", assembly)
    elseif isnothing(format) || format == "mmCIF"
        _wget(pdb_id, selection_function; format="cif", assembly)
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

    # biological assembly: 3CNA's asymmetric unit is a single chain, but its
    # assembly 1 is the crystallographic tetramer (see richards.jl SASA test)
    tetramer = wget("3CNA", "protein"; assembly=1)
    @test Set(unique(chain.(tetramer))) == Set(["A", "A-2", "A-3", "A-4"])
    @test length(tetramer) == 7228
    tetramer_pdb = wget("3CNA", "protein"; format="PDB", assembly=1)
    @test length(tetramer_pdb) == length(tetramer)
    @test_throws "identifier code or internet connection" wget("####") 
end
