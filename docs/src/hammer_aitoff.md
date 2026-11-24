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

Here, we use the capsid of the Zika virus, in a coarse-grained (SIRAH) representation.
The capsid is formed by 60 copies of an asymmetric unit (segments `S1` to `S60`), each
with three copies of the envelope (E) protein (chains `K`, `M` and `O`) and three copies 
of the membrane (M) protein (chains `L`, `N` and `P`). We keep only the `GC` beads, which 
are at the positions of the ``C_\alpha`` atoms, such that there is one bead per residue: 
103,680 beads.

[`hammer_aitoff_projection`](@ref) projects the atoms, and returns a [`HammerAitoffProjection`](@ref),
which contains the atoms, the center of the projection and the projected coordinates of each
atom (fields `x` and `y`). The input structure is not modified.

```@example hammer
using PDBTools, LazyArtifacts
capsid = read_pdb(PDBTools.test_file("capsid.pdb"), "name GC")
proj = hammer_aitoff_projection(capsid)
```

```@example hammer
proj.x[1:3], proj.y[1:3]
```

## Heat maps of per-atom properties

[`hammer_aitoff_map`](@ref) projects the atoms, divides the map into `nbins × nbins` bins and
computes, for each bin, a statistic of the values of the atoms that fell inside it. Any per-atom
property can be used: B-factors, partial charges, contributions to coordination numbers,
distances, etc. The resulting [`HammerAitoffMap`](@ref) is plotted with `heatmap` (a Makie
backend must be loaded, e.g. `using CairoMakie`).

Here, each bead is colored by its distance to the center of the capsid, which shows the
arrangement of the E proteins on the surface of the virus, repeated by the icosahedral 
symmetry of the capsid:

```@example hammer
using CairoMakie
c = center_of_mass(capsid)
dist = [ sqrt((at.x - c[1])^2 + (at.y - c[2])^2 + (at.z - c[3])^2) for at in capsid ]
map = hammer_aitoff_map(capsid, dist; nbins = 100)
```

```@example hammer
heatmap(map; colorbar_title = "distance to center (Å)")
```

Bins without atoms (for example, outside the ellipse, or at the small openings of the capsid
at its five-fold symmetry axes) have no value (`NaN`), and are not colored.

The `statistic` keyword of `hammer_aitoff_map` defines how the values inside a bin are combined:

- `:mean` (default): for intensive properties, such as B-factors or distances;
- `:sum`: for additive properties, such as atomic contributions to a coordination number;
- `:count`: number of atoms in each bin (the values are ignored).

To compare several structures on the same color scale, fix the limits with the `clims` 
keyword of `heatmap`, e.g. `clims = (170, 240)`.

The bin centers and the binned values are available in the `xbin`, `ybin` and `values` fields
of the map, where `values[j, i]` is the value of the bin with center `(xbin[i], ybin[j])`:

```@example hammer
size(map.values)
```

A map can also be computed from a projection, such that maps of several properties of the
same atoms do not require projecting the atoms again. Here, the number of beads in each
bin, which depends on how much protein is found along each direction from the center (it is
larger, on average, where the M proteins lie beneath the E proteins):

```@example hammer
map_count = hammer_aitoff_map(proj, ones(length(capsid)); nbins = 100, statistic = :count)
heatmap(map_count; colorbar_title = "number of beads")
```

## Contours of chains and groups of atoms

With `contours = :chain`, a line is drawn around the region occupied by each chain. With 
`contours = :group`, the atoms are grouped by a function that receives each atom and returns
the name of its group. Here, the atoms are grouped by segment, which shows the 60 asymmetric 
units of the capsid:

```@example hammer
heatmap(map; 
    contours = :group, groupfun = segname, contour_color = :white, 
    colorbar_title = "distance to center (Å)",
)
```

Any property of the atoms can define the groups. For example, the chains of a structure with 
many copies of a few proteins can be grouped by protein:

```julia
groups = Dict("A" => "protein 1", "B" => "protein 1", "C" => "protein 2")
heatmap(map; contours = :group, groupfun = at -> groups[chain(at)])
```

The figure can be saved with `CairoMakie.save("map.png", fig)`.

## Reference

```@docs
hammer_aitoff_projection
HammerAitoffProjection
hammer_aitoff_map
HammerAitoffMap
heatmap(::HammerAitoffMap)
```
