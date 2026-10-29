using TestItemRunner

@testmodule AllocTest begin
    # This module defines the Allocs struct and the comparison operators
    # to conditionally compare the number of allocations based on the
    # BUILD_IS_PRODUCTION_BUILD environment variable.
    export Allocs
    @kwdef struct Allocs
        prodbuild::Bool = haskey(ENV, "BUILD_IS_PRODUCTION_BUILD") && ENV["BUILD_IS_PRODUCTION_BUILD"] == "true"
        allocs::Int
    end
    Allocs(allocs::Int) = Allocs(; allocs)
    import Base: ==, >, <
    ==(a::Real, b::Allocs) = b.prodbuild ? Int(a) == b.allocs : true
    <(a::Real, b::Allocs) = b.prodbuild ? Int(a) < b.allocs : true
    ==(a::Allocs, b::Real) = a.prodbuild ? a.allocs == Int(b) : true
    <(a::Allocs, b::Real) = a.prodbuild ? a.allocs < Int(b) : true
end

@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    using PDBTools: PDBToolsBase
    # PDBTools adds methods to the functions and types of PDBToolsBase (e.g. for selection
    # strings, and for the secondary structure functions with vectors of atoms)
    Aqua.test_all(PDBTools;
        piracies=(treat_as_own=[
            PDBToolsBase.Atom, PDBToolsBase.read_pdb, PDBToolsBase.read_mmcif, PDBToolsBase.write_pdb,
            PDBToolsBase.write_mmcif, PDBToolsBase.read_xyz, PDBToolsBase.write_xyz,
            PDBToolsBase.positions, PDBToolsBase.maxmin,
        ],),
    )
end
