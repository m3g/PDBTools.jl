using Documenter
using Plots
using PDBTools
using LazyArtifacts # loads the extension that defines PDBTools.test_file
# Download the test files before the doctests, so that the download messages are not part of their output
PDBTools.test_file("structure.pdb") # files of PDBToolsBase
PDBTools.test_file("hbonds.pdb") # files of PDBTools
ENV["GKSwstype"] = "100"
ENV["LINES"] = 10
ENV["COLUMNS"] = 120
ENV["GKSwstype"] = "nul"
makedocs(
    modules = [
        PDBTools, 
        PDBTools.PDBToolsBase,
        isdefined(Base, :get_extension) ? Base.get_extension(PDBTools, :Plotting) : PDBTools.Plotting
    ],
    # The core functions are defined in PDBToolsBase: the docstrings of the public names of 
    # both modules must be included in the manual.
    checkdocs = :public,
    # Source links of the docstrings of PDBToolsBase
    remotes = Dict(pkgdir(PDBTools.PDBToolsBase) => (Documenter.Remotes.GitHub("m3g", "PDBToolsBase.jl"), "main")),
    # The visualization page embeds the structures in the HTML
    format = Documenter.HTML(top_menu = true, size_threshold_ignore = ["visualization.md"]),
    sitename = "PDBTools.jl",
    pages = [
        "Home" => "index.md",
        "Migrating to 4.0" => "migration.md",
        "Getting started" => Any[
            "Read and Write" => "readwrite.md",
            "Structures and unit cells" => "structure.md",
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
            "Coordinate manipulations" => "coordinates.md",
            "Visualization" => "visualization.md",
        ],
        "Transfer Free Energy" => Any[
            "Transfer Free Energy" => "mvalues/mvalue.md",
            "Protein Denaturation" => "mvalues/protein_denaturation.md",
            "Conformational Changes" => "mvalues/conformational_changes.md",
            "Cofactor Binding" => "mvalues/cofactor_binding.md",
            "Accessibility Model" => "mvalues/accessibility.md",
            "Record Model" => "mvalues/record.md",
            "SASA precomputation" => "mvalues/sasa_precomputing.md",
            "Alternative Calculations" => "mvalues/alternative_mvalues.md",
        ],
        "Other" => Any[
            "Auxiliary functions" => "auxiliary.md",
            "Examples" => "examples.md",
        ],
    ],
)
deploydocs(
    repo = "github.com/m3g/PDBTools.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
