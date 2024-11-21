import Pkg
Pkg.add("Documenter")
using Documenter
using PDBTools
makedocs(
    modules = [PDBTools],
    sitename = "PDBTools.jl",
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Read and Write" => "readwrite.md",
        "Selections" => "selections.md",
        "Element properties" => "elements.md",
        "Auxiliary functions" => "auxiliary.md",
        "Examples" => "examples.md",
        "Tools" => "tools.md",
        "Help entries" => "help.md",
    ],
)
deploydocs(
    repo = "github.com/m3g/PDBTools.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
