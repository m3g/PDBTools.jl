import Pkg
Pkg.add("Documenter")
using Documenter
using Plots
using PDBTools
ENV["LINES"] = 10
ENV["COLUMNS"] = 120
makedocs(
    modules = [
        PDBTools, 
        isdefined(Base, :get_extension) ? Base.get_extension(PDBTools, :Plotting) : PDBTools.Plotting
    ],
    sitename = "PDBTools.jl",
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Read and Write" => "readwrite.md",
        "Selections" => "selections.md",
        "Iterators" => "iterators.md",
        "Element properties" => "elements.md",
        "Contact maps" => "contacts.md",
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
