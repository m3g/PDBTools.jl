```@meta
CollapsedDocStrings = true
```

# [Visualization](@id visualization)

The `visualize` function returns an interactive 3D view of a structure, rendered with 
[3Dmol.js](https://3dmol.csb.pitt.edu). The view is displayed in environments that render HTML: the 
VS Code plot pane, Pluto and Jupyter notebooks, and Documenter pages (like this one). In the Julia REPL, the view is opened in the default web browser. 
From scripts, use `open_browser(view)`, or save the view to an HTML file with `save("view.html", view)`.

The secondary structure shown in the cartoon representation is computed by 3Dmol.js. The 3Dmol.js 
library is loaded from a CDN, so an internet connection is required to display the views.

```@docs
visualize
PDBTools.StructureView
save(::AbstractString, ::PDBTools.StructureView)
open_browser
```

## Examples

By default, proteins and nucleic acids are shown as cartoons, other molecules (ligands, ions) as sticks, 
and water molecules are hidden (unless there is nothing else to show, or `water=true`). Hovering over
atoms shows their residue, chain, and name:

```@example visualization
using PDBTools
atoms = read_pdb(PDBTools.TESTPDB, "protein")
visualize(atoms)
```

Secondary structure coloring, with a transparent molecular surface:

```@example visualization
visualize(atoms; color=:ss, surface=true)
```

Atomistic representations of a selection:

```@example visualization
visualize(atoms, "residue <= 20"; style=:ballandstick)
```

Coloring the atoms by any property, with one value per atom. Here, the solvent accessible
surface area of each atom:

```@example visualization
s = sasa_particles(atoms)
visualize(atoms; color_by=[s[i] for i in eachindex(atoms)], colormap=:roygb)
```

The edges of the periodic box can be drawn by providing the unit cell, as a 3x3 matrix with the lattice
vectors as columns (as returned by `read_unitcell`) or as a vector of box sides. By default the box
origin is at `(0, 0, 0)`; use `unitcell_origin=:center` for boxes centered on the system:

```julia
atoms = read_pdb(PDBTools.TESTPBC)
visualize(atoms; unitcell=read_unitcell(PDBTools.TESTPBC))
visualize(atoms; unitcell=[85.0, 85.0, 85.0], unitcell_origin=:center)
```
