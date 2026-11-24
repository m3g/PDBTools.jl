module TestFilesExt

using LazyArtifacts: LazyArtifacts, @artifact_str
import PDBTools: test_file
import PDBTools

function test_file(name::AbstractString)
    name in PDBTools._extra_test_files && return joinpath(artifact"pdbtools_extra_test_files", name)
    name in PDBTools._zika_test_files && return joinpath(artifact"zika_test_files", name)
    # The pdbtools set (also used by PDBToolsBase), from the Artifacts.toml of PDBTools
    file = joinpath(artifact"pdbtools_test_files", name)
    isfile(file) || throw(ArgumentError("Test file \"$name\" not found."))
    return file
end

end
