```@meta
CollapsedDocStrings = true
```

# Hammer-Aitoff projection

The Hammer-Aitoff projection maps the surface of a sphere onto an ellipse, preserving areas. Applied to a molecular structure, it gives a two-dimensional map of the whole surface of the system, as seen from its center, in a single figure. This is particularly useful for large, roughly spherical assemblies, such as viral capsids and envelopes, where any single view of the 3D structure hides half of the surface.

## How the projection works

For each atom, the direction from a center (by default, the center of mass of the selection) to the atom is converted into a latitude ``\varphi \in [-\pi/2, \pi/2]`` and a longitude ``\lambda \in [-\pi, \pi]``:

```math
\varphi = \arcsin\left(\frac{z}{r}\right), \qquad \lambda = \operatorname{atan}(y, x)
```

where `(x, y, z)` are the coordinates of the atom relative to the center and `r` its distance to the center. These angles are then mapped onto the plane with

```math
x_H = \frac{2\sqrt{2}\,\cos\varphi\,\sin(\lambda/2)}{\sqrt{1 + \cos\varphi\,\cos(\lambda/2)}}, \qquad
y_H = \frac{\sqrt{2}\,\sin\varphi}{\sqrt{1 + \cos\varphi\,\cos(\lambda/2)}}
```

so that all points lie inside the ellipse ``(x_H / 2\sqrt{2})^2 + (y_H / \sqrt{2})^2 \leq 1``. Because the projection preserves areas, regions of the same size on the surface occupy the same area in the map, regardless of where they are.

Only the direction of each atom is used: atoms along the same line from the center fall on the same point of the map. The value shown in each region of the map is a statistic (mean, sum or count) of the values of all atoms projected onto it.

## Projecting the atoms

[`hammer_aitoff_projection`](@ref) returns the projected coordinates of each atom. The input structure is not modified.

```@example hammer
using PDBTools, LazyArtifacts
atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
x, y = hammer_aitoff_projection(atoms)
x[1:3], y[1:3]
```

## Heat maps of per-atom properties

[`hammer_aitoff_map`](@ref) projects the atoms, divides the map into `nbins × nbins` bins and computes, for each bin, a statistic of the values of the atoms that fell inside it. Any per-atom property can be used: B-factors, partial charges, contributions to coordination numbers, distances, etc. The resulting [`HammerAitoffMap`](@ref) is plotted with `heatmap` (`Plots` must be loaded).

Here, each atom is colored by its distance to the center of mass, which highlights protrusions and cavities of the surface:

```@example hammer
using Plots
c = center_of_mass(atoms)
dist = [sqrt((at.x - c[1])^2 + (at.y - c[2])^2 + (at.z - c[3])^2) for at in atoms]
map = hammer_aitoff_map(atoms, dist; nbins = 100)
```

```@example hammer
heatmap(map; colorbar_title = "distance to center (Å)")
```

The `statistic` keyword of `hammer_aitoff_map` defines how the values inside a bin are combined:

- `:mean` (default): for intensive properties, such as B-factors or distances;
- `:sum`: for additive properties, such as atomic contributions to a coordination number;
- `:count`: number of atoms in each bin (the values are ignored).

To compare several structures on the same color scale, fix the limits with `clims`:

```@example hammer
heatmap(map; clims = (0, 40), colorbar_title = "distance to center (Å)")
```

The bin centers and the binned values are available in the `xbin`, `ybin` and `values` fields of the map,
where `values[j, i]` is the value of the bin with center `(xbin[i], ybin[j])`:

```@example hammer
size(map.values)
```

## Contours of chains and groups of chains

With `contours = :chain`, a line is drawn around the region occupied by each chain:

```@example hammer
heatmap(map; contours = :chain, colorbar_title = "distance to center (Å)")
```

With `contours = :group`, chains are grouped by a function that receives the chain name and returns the name of its group. This is useful when the structure has many copies of a few proteins (for example, the capsid and envelope proteins of a virus):

```julia
groups = Dict("A" => "protein 1", "B" => "protein 1", "C" => "protein 2")
heatmap(map; contours = :group, groupfun = ch -> groups[ch])
```

The figure can be saved with `savefig`.

## Reference

```@docs
hammer_aitoff_projection
hammer_aitoff_map
HammerAitoffMap
heatmap(::HammerAitoffMap)
```
