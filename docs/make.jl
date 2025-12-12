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
        "Read and Write" => "readwrite.md",
        "Selections" => "selections.md",
        "Iterators" => "iterators.md",
        "Atom and element properties" => "elements.md",
        "Contact maps" => "contacts.md",
        "Dihedrals and Ramachandran" => "dihedrals.md",
        "Secondary structures" => "secondary_structure.md",
        "Hydrogen bonds" => "hydrogen_bonds.md",
        "Solvent Accessible Area" => "sasa.md",
        "Transfer Free Energy" => "mvalue.md",
        "Auxiliary functions" => "auxiliary.md",
        "Examples" => "examples.md",
        "Tools" => "tools.md",
        "Help entries" => "help.md",
        "Experimental" => Any[ 
            "Structure" => "structure.md",
        ]
    ],
)
deploydocs(
    repo = "github.com/m3g/PDBTools.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
