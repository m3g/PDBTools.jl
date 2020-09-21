import Pkg
Pkg.add("Documenter")
using Documenter
using PDBTools
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[PDBTools],
    sitename="PDBTools.jl",
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Read and Write" => "readwrite.md",
        "Selections" => "selections.md",
        "Element properties" => "elements.md",
        "Auxiliary functions" => "auxiliary.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/PDBTools.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)
