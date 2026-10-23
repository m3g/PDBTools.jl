module TestFilesExt

using LazyArtifacts: LazyArtifacts, @artifact_str
import PDBTools: test_file
import PDBTools

function test_file(name::AbstractString)
    # The files shared with PDBToolsBase are provided by PDBToolsBase.test_file
    name in PDBTools._extra_test_files && return joinpath(artifact"pdbtools_extra_test_files", name)
    return PDBTools.PDBToolsBase.test_file(name)
end

end
