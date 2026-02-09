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
    Aqua.test_all(PDBTools)
end

#@testitem "Doctests" begin
#    using Documenter: doctest
#    using PDBTools
#    doctest(PDBTools)
#end