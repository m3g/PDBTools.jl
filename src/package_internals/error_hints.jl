#
# Error hints for functions whose methods are defined in package extensions: if the
# package that triggers the extension was not loaded, a MethodError is thrown, and
# the hint tells the user which package to load.
#
const _extension_functions = (
    (f = plot_ss_diagram, ext = :MakiePlotting, requires = "a Makie backend", pkg = "CairoMakie"),
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
    exc = MethodError(PDBTools.plot_ss_diagram, (1,))
    msg = sprint(PDBTools._extension_error_hint, exc, (Int,), (;))
    if isnothing(Base.get_extension(PDBTools, :MakiePlotting))
        @test occursin("using CairoMakie", msg)
    else
        @test isempty(msg)
    end
    # unrelated MethodErrors get no hint
    @test isempty(sprint(PDBTools._extension_error_hint, MethodError(sin, ("a",)), (String,), (;)))
end
