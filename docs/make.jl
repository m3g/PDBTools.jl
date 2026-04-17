using Documenter
using Plots
using PDBTools
ENV["GKSwstype"] = "100"
ENV["LINES"] = 10
ENV["COLUMNS"] = 120
makedocs(
    modules = [
        PDBTools, 
        isdefined(Base, :get_extension) ? Base.get_extension(PDBTools, :Plotting) : PDBTools.Plotting
    ],
    format = Documenter.HTML(top_menu = true),
    sitename = "PDBTools.jl",
    pages = [
        "Home" => "index.md",
        "Getting started" => Any[
            "Read and Write" => "readwrite.md",
            "Selections" => "selections.md",
            "Iterators" => "iterators.md",
            "Atomic and molecular properties" => "elements.md",
        ],
        "Structure" => Any[
            "Contact maps" => "contacts.md",
            "Dihedrals and Ramachandran" => "dihedrals.md",
            "Secondary structures" => "secondary_structure.md",
            "Hydrogen bonds" => "hydrogen_bonds.md",
            "Solvent Accessible Area" => "sasa.md",
            "Coordinate manipulations" => "coordinates.md"
        ],
        "m-values" => Any[
            "Transfer Free Energy" => "mvalue.md",
        ],
        "Other" => Any[
            "Auxiliary functions" => "auxiliary.md",
            "Examples" => "examples.md",
            "Experimental" => Any[ 
                "Structure" => "structure.md",
            ],
        ],
    ],
)
deploydocs(
    repo = "github.com/m3g/PDBTools.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
