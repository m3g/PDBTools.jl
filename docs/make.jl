#using DocumenterLaTeX
using Documenter
using PDBTools
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[PDBTools],
    sitename="PDBTools.jl",
    pages = [
        "Installation" => "index.md",
        "Read and Write" => "readwrite.md",
        "Selections" => "selections.md",
        "Auxiliary functions" => "auxiliary.md",
    ]
)



