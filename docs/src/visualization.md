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

## Groups with different representations

Different groups of atoms can be shown with different representations in the same view. Each group is
given by a pair of the atoms and a `NamedTuple` with the representation options (`style`, `color`, 
`color_by`, `colormap`, `color_range`, `ligands`, `water`, `surface`, `opacity`, and `selection`). The groups
can be selections of the same vector of atoms. Here, a cartoon of the protein, and the
acidic and basic residues as sticks:

```@example visualization
visualize(atoms, 
    "protein" => (color="white",), 
    "acidic" => (style=:sticks, color="red"),
    "basic" => (style=:sticks, color="blue"),
)
```

Or the groups can be different vectors of atoms. Each group is a separate model in the view, such that
bonds are not computed between atoms of different groups. For example, a set of points around the
protein can be shown as dots (`style=:dots`), colored by a property: 

```@example visualization
points = [ 
    Atom(name="X", resname="PNT", resnum=i, x=at.x + 3 * cos(i), y=at.y + 3 * sin(i), z=at.z, beta=sin(i)^2) 
    for (i, at) in enumerate(select(atoms, "name CA")) 
]
visualize(
    atoms => (color="white",),
    points => (style=:dots, color_by=beta.(points), colormap=:rwb, color_range=(1, 0)),
)
```

The `opacity` option sets the opacity of the atoms of each group, from `0` (invisible) to `1` (opaque). Here,
the protein is shown as a transparent cartoon:

```@example visualization
visualize(atoms, "protein" => (color="white", opacity=0.4), "acidic" => (style=:ballandstick, color="red"))
```

The options that apply to the whole view (`hover`, `unitcell`, `unitcell_origin`, `unitcell_color`,
`background`, `width`, and `height`) are given as keyword arguments.

## Periodic boxes

The edges of the periodic box can be drawn by providing the unit cell, as a 3x3 matrix with the lattice
vectors as columns (as returned by `read_unitcell`) or as a vector of box sides. By default the box
origin is at `(0, 0, 0)`; use `unitcell_origin=:center` for boxes centered on the system:

```julia
atoms = read_pdb(PDBTools.TESTPBC)
visualize(atoms; unitcell=read_unitcell(PDBTools.TESTPBC))
visualize(atoms; unitcell=[85.0, 85.0, 85.0], unitcell_origin=:center)
```
