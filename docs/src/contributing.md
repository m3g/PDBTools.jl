# Contributing

Contributions are welcome! This page summarizes the conventions that new functions should 
follow, such that they work well with the rest of the package. Following them makes the 
review of a pull request much faster.

## Interface

### Inputs are atoms

Functions that operate on structures receive the atoms, as an `AbstractVector{<:Atom}`. This
accepts the result of [`read_pdb`](@ref) (a [`Structure`](@ref)), selections, residues, 
chains, and any other vector of atoms. Do not accept file names or selection strings: reading 
and selecting are done before, with the functions of the package:

```julia
atoms = read_pdb("file.pdb", "protein and chain A")   # reading and selecting
result = new_function(atoms)                           # the new function receives atoms
```

Name the argument `atoms`, as in the rest of the package.

### Reuse the functions of the package

Before writing a helper function, check if the package already provides it. For example:

- per-residue values: residues, chains, segments, models and structures are vectors of atoms,
  so `sum(f, residue)` (or `mean(f, residue)`, with `Statistics`) computes the value of a 
  residue from any atomic property, and [`eachresidue`](@ref) iterates over the residues;
- sequences: [`getseq`](@ref), [`oneletter`](@ref);
- geometry: [`center_of_mass`](@ref), [`positions`](@ref), [`distance`](@ref);
- solvent accessible surface areas: [`sasa_particles`](@ref) and [`sasa`](@ref).

### Compute first, then plot

Computation and plotting are separate. A function computes the result and returns it as a
custom type, which carries everything needed to plot it, and can be used without any 
plotting package. The plot is a method, for that type, of a function of the plotting package
(`heatmap`, `plot`, `scatter` of Makie):

| Computation | Result | Plot |
|:------------|:-------|:-----|
| [`contact_map`](@ref) | [`ContactMap`](@ref) | `heatmap(map)` |
| [`Ramachandran`](@ref) | `Ramachandran` | `scatter(ram)` |
| [`hammer_aitoff_map`](@ref) | [`HammerAitoffMap`](@ref) | `heatmap(map)` |
| [`helical_wheel`](@ref) | [`HelicalWheel`](@ref) | `plot(wheel)` |
| [`secondary_structure_diagram`](@ref) | [`SecondaryStructureDiagram`](@ref) | `plot(diagram)` |

The result types have a `show` method, and fields with descriptive names.

The plots are made with Makie, in the `MakiePlotting` extension (in `ext/MakiePlotting`), which
is loaded with any Makie backend (e.g. `using CairoMakie`). The plotting methods are methods of
the functions of Makie (`plot`, `heatmap`, `scatter`), for the type returned by the computation.
Each type has three methods (for `f` = `plot`, `heatmap` or `scatter`):

- `f(x; kwargs...)`: the complete figure, with its color bar, returned as a `Makie.Figure`.
  This is the method described in the section of the manual of the feature. It has the 
  `size`, `font` and `fontsize` keywords, and the keywords of the position method.
- `f(fig[i, j], x; kwargs...)`: draws the plot and its color bar into a position of a figure,
  and returns the `GridLayout` of the position. 
- `f!(ax, x; kwargs...)`: draws the plot into an existing `Axis`, without color bar and 
  without changing the axis, and returns the plot object.

The last two methods are documented only in the [Advanced plotting](@ref) section of the manual.

Plotting functions return the figure, and do not save it: users save figures with 
`CairoMakie.save` (so there are no `output` or `outfile` keywords). The plots with `Plots.jl`
of the `PlotsPlotting` extension are deprecated, and do not receive new features.

### One quantity per function

If two quantities are computed, prefer two functions, each returning one quantity, to a 
function returning a tuple. For example, [`hydrophobic_moment`](@ref) returns the moment (a 
vector) and [`mean_hydrophobicity`](@ref) returns a number.

### Periodic boundary conditions

Functions that support periodic boundary conditions have a `unitcell` keyword argument, which
is `false` (default, no periodic boundary conditions), `true` (use the unit cell of the 
`Structure`), or a [`UnitCell`](@ref). See [Structures and unit cells](@ref).

## Names

- Function names are descriptive, in `snake_case`, and avoid abbreviations: 
  `hammer_aitoff_projection`, not `hammer`; `secondary_structure_diagram`, not `ss`. Names 
  are part of the public interface, and are hard to change after a release.
- Type names are in `CamelCase` (`HammerAitoffMap`).
- Export only the functions and types that users call. Helper functions are internal, and 
  their names start with `_` (e.g. `_binned_statistic_2d`). Constants and data tables (e.g. 
  hydrophobicity scales) are internal, and their sources are cited in the documentation.

## Code organization

- Source files go in `src` (or in a subdirectory of `src`), with the `export` statements at the
  top of the file. They are included in `src/PDBTools.jl`, in the section they belong to, 
  before the precompilation workload.
- Plotting code goes in the `ext/MakiePlotting` directory, in a `<feature>.jl` file included
  in `ext/MakiePlotting/MakiePlotting.jl`.
- Dependencies: avoid new dependencies. A heavy dependency, needed by a single feature, should
  be a weak dependency, with a package extension. New dependencies and extensions must be 
  clearly stated in the description of the pull request.

## Tests

- Tests are `@testitem`s, in the same file as the code they test (plotting tests in the files
  of the `ext` directory). They are run with `]test`, or individually in VS Code.
- Test the results (with known values, or with properties that must hold), and also the 
  errors, with `@test_throws`, for invalid inputs.
- Test files are obtained with `PDBTools.test_file("name")`, after `using LazyArtifacts`. The
  test files are lazy artifacts, from the [MolSimTestFiles](https://github.com/m3g/MolSimTestFiles)
  repository, and are not part of the package. If a new test file is needed, mention it in 
  the pull request: it will be added to a new release of MolSimTestFiles.
- The quality tests of [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) must pass.

## Documentation

- Every exported function and type has a docstring, with the signature, the description of 
  the arguments and keyword arguments, the return value, and an example. The examples start 
  with `using PDBTools, LazyArtifacts`, and use `PDBTools.test_file`. The documentation build
  checks that every public name is in the manual.
- New features are described in a page of `docs/src`, which must be added to the `pages` of 
  `docs/make.jl` (otherwise it is not reachable in the documentation). The page explains what
  the feature is for, and contains the docstrings in a `@docs` block.
- Examples in the pages are `@example` blocks, which are run when the documentation is built,
  such that the results and figures are shown and the examples are known to work. Use 
  ```` ```julia ```` blocks only for code that cannot be run (e.g. code that needs data that 
  is not available).
- Equations are written in ```` ```math ```` blocks, or inline as ``` ``\alpha`` ```.
- The documentation is built with `julia --project=docs docs/make.jl`. For pull requests from 
  branches of this repository, a preview is deployed to 
  `https://m3g.github.io/PDBTools.jl/previews/PR<number>/`. For pull requests from forks, 
  the documentation is built (and checked), but not deployed: build it locally to see it.

## Changelog

Every pull request adds an entry to `CHANGELOG.md` (this is checked in the CI), in the 
section of the version under development (at the top), below its header line. The entries 
start with a badge, which classifies the change:

```
- ![FEATURE][badge-feature] `new_function` computes ...
```

The badges are `FEATURE`, `ENHANCEMENT`, `BUGFIX`, `BREAKING`, `DEPRECATION`, `EXPERIMENTAL`
and `INFO`. The entry describes the change for users: what can be done, with which functions.

New features are released in a new minor version (e.g. 4.3.0), and bug fixes in a new patch 
version (e.g. 4.2.1). Breaking changes require a new major version, and are avoided.
