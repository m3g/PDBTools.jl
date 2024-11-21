using TestItemRunner
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