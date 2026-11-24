```@meta
CollapsedDocStrings = true
```
# Contact and distance maps

The `contact_map` function computes a contact map for a structure or a pair of structures.
These structures are typically proteins, but any structures defined by sequences of residues
can be provided as inputs. The maps are plotted with `heatmap` (a Makie backend must be
loaded, e.g. `using CairoMakie`).

```@docs
contact_map
heatmap(::ContactMap)
```

!!! note
    The distance computed in these functions is the **minimum distance** between the 
    atoms of the residues. If, for example, a contact map computed from the distance
    between Cα atoms is desired, select the atoms before computing the map:
    ```julia
    cA = select("chain A and name CA")
    contact_map(cA; dmax=8.0)
    ```
    Importantly, the maximum distance defining a contact has to be adjusted, as
    it is, by default, `dmax=4.0`, which is a reasonable contact measure from
    minimum distance between atoms, but too short for distances between backbone
    atoms.

## Contact map

A typical usage consists in computing the contact map and plotting it:

```@example contacts
using PDBTools, LazyArtifacts, CairoMakie
ats = read_pdb(PDBTools.test_file("dimer.pdb"));
cA = select(ats, "chain A");
cB = select(ats, "chain B");
map = contact_map(cA, cB) # contact map between chains A and B
heatmap(map)
```

## Distance map

In the example above we opted to plot a discrete contact map, with the default contact 
distance `dmax=4.0`. Now we change two parameters: `discrete=false` and `dmax=12.0`, to
compute a distance map up to a greater distance:

```@example contacts
map = contact_map(cA, cB; discrete=false, dmax=8.0) # contact map between chains A and B
heatmap(map)
```

## Single structure

Similarly, we can produce plots for the contact map of a single structure. Here, we 
showcase the use of the `gap` parameter, to ignore residues closer in the sequence
by less than 4 residues:

```@example contacts
distance_map = contact_map(cA; gap=4, discrete=false, dmax=8.0) # chain A only
heatmap(distance_map; colormap=:davos)
```

```@example contacts
discrete_map = contact_map(cA; gap=4, discrete=true, dmax=8.0) # chain A only
heatmap(discrete_map)
```

## Difference (and sum) of maps

Contact maps of the same type (discrete *or* continuous), obtained for the same sequence, but with
different conformations, can be compared by subtraction or summation. 

For example, here we load two models from a PDB file that contains multiple conformations of a 
protein, in solution, and compute the contact maps of the two models, and plot the difference
of the maps. 
```@example contacts
pdb = wget("2cpb", "model 1 2")
models = collect(eachmodel(pdb))
c1 = contact_map(models[1])
c2 = contact_map(models[2])
heatmap(c2 - c1)
```
Blue dots indicate contacts present in the first set, but not in the second, and red dots
the contacts present in the second but not in the first.

A difference of distance maps can be similarly obtained by computing continuous contact maps:
```@example contacts
c1 = contact_map(models[1]; discrete=false, dmax=12.0)
c2 = contact_map(models[2]; discrete=false, dmax=12.0)
c_diff = c2 - c1
heatmap(c_diff)
```

!!! note
    Contacts farther than the tolerance set are `missing`. When summing or subtracting 
    distance maps, `missing` values are propagated. For discrete maps, `true` contacts
    are preserved, resulting in `+1` or `-1` depending on the contact being present in 
    the first or second structure.

## Customizing the plot

Common customization options of `heatmap` are (see [`heatmap(::ContactMap)`](@ref) for all of them):

- `xstep`, `ystep`: the stride of the residue labels on the axes;
- `colormap`: the color map (default: `:grayC` for distances, black for binary maps);
- `clims`: the range of the color scale;
- `colorbar`, `colorbar_title`: the color bar and its title (default: "distance (Å)" for
  distances, no color bar for binary maps);
- `size`: the size of the figure.

Above, we illustrated this by changing the color map of the distance map. The figure can be
saved with `CairoMakie.save("map.png", fig)`.

## REPL visualization of the contact matrix

Although this is not considered a fundamental feature, a quick visualization of the 
contact matrix can be obtained in the REPL by showing the sparse matrix representation
directly:

```@example contacts
c1.matrix
```

## Indexing

The `ContactMap` data structure can be indexed to extract the contacts of a specific 
residue. For example:

```@example getindex_map
using PDBTools, LazyArtifacts
ats = read_pdb(PDBTools.test_file("dimer.pdb"));
cA = select(ats, "chain A");
cB = select(ats, "chain B");
map = contact_map(cA, cB; discrete=false, dmax=12.0)
map[235,7] # Distance of 235 to 7
```

Or slicing:
```@example getindex_map
map[235,:] # all distances below 12.0 Angs of residue 235 of cA with cB
```

## Data structure and auxiliary functions

```@docs
ContactMap
```
