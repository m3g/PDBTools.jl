# Deprecated: plots with Plots.jl

!!! warning
    The plots of PDBTools are made with Makie (see, for example, [Contact and distance maps](@ref) and
    [Dihedrals and Ramachandran plots](@ref)). The plots with `Plots.jl` described here are 
    deprecated: they are kept for compatibility, without new features, and will be removed
    in a future breaking release.

Previous versions of PDBTools provided plots of contact maps and of Ramachandran plots with
`Plots.jl`. These are still available when `Plots` is loaded (and no Makie backend is
used in the same session, because both packages export `heatmap` and `scatter`):

```julia
using PDBTools, Plots
map = contact_map(atoms)
heatmap(map) # Plots.heatmap(::ContactMap)

ram = Ramachandran(atoms)
scatter(ram) # Plots.scatter(::Ramachandran)
```

The keywords of these functions are those of `Plots.heatmap` and `Plots.scatter`, with the 
following additional options for contact maps: `xstep` and `ystep` (stride of the residue
labels on the axes), `color` (the color map), `clims` (the range of the color scale) and
`colorbar_title`. The figures are saved with `savefig`.

To update code that uses these functions, replace `using Plots` by `using CairoMakie`, the
keyword `color` by `colormap`, and `savefig(plt, file)` by `CairoMakie.save(file, fig)`.
