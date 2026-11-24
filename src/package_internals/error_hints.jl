#
# Error hints for functions whose methods are defined in package extensions: if the
# package that triggers the extension was not loaded, a MethodError is thrown, and
# the hint tells the user which package to load.
#
#
# (The plotting functions are methods of functions of the plotting packages, as
# `Makie.plot`, so they cannot have hints here.)
#
const _extension_functions = (
    (f = test_file, ext = :TestFilesExt, requires = "the LazyArtifacts package", pkg = "LazyArtifacts"),
)

function _extension_error_hint(io, exc, argtypes, kwargs)
    for (; f, ext, requires, pkg) in _extension_functions
        if exc.f === f && isnothing(Base.get_extension(@__MODULE__, ext))
            print(io, "\n\n`PDBTools.$(nameof(f))` requires ")
            printstyled(io, requires; color = :cyan, bold = true)
            print(io, ". Load it with, for example:\n\n    using $pkg\n")
        end
    end
end

function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(_extension_error_hint, MethodError)
    end
end

@testitem "extension error hints" begin
    using PDBTools
    exc = MethodError(PDBTools.test_file, ("structure.pdb",))
    msg = sprint(PDBTools._extension_error_hint, exc, (String,), (;))
    if isnothing(Base.get_extension(PDBTools, :TestFilesExt))
        @test occursin("using LazyArtifacts", msg)
    else
        @test isempty(msg)
    end
    # unrelated MethodErrors get no hint
    @test isempty(sprint(PDBTools._extension_error_hint, MethodError(sin, ("a",)), (String,), (;)))
end
