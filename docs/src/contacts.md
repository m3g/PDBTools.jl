```@meta
CollapsedDocStrings = true
```
# Contact and distance maps

The `contact_map` function computes a contact map for a structure or a pair of structures.
These structures are typically proteins, but any structures defined by sequences of residues
can be provided as inputs. The `heatmap` function, from `Plots`, is overloaded here
to provide a convenient way to plot the contact or distance maps.

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
using PDBTools
using Plots
ats = read_pdb(PDBTools.DIMERPDB);
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
map = contact_map(cA, cB; discrete=false, dmax=12.0) # contact map between chains A and B
heatmap(map)
```

## Single structure

Similarly, we can produce plots for the contact map of a single structure. Here, we 
showcase the use of the `gap` parameter, to ignore residues closer in the sequence
by less than 4 residues:

```@example contacts
distance_map = contact_map(cA; gap=4, discrete=false, dmax=12.0) # chain A only
discrete_map = contact_map(cA; gap=4, discrete=true, dmax=12.0) # chain A only
plot(
    heatmap(distance_map; colorbar=nothing, color=:davos), 
    heatmap(discrete_map); 
    layout=(1,2), size=(800,500)
)
```

## Difference (and sum) of maps

Contact maps of the same type (discrete *or* continous), obtained for the same sequence, but with
different conformations, can be compared by subtraction or summation. 

For example, here we load two models from a PDB file that contains multiple conformations of a 
protein, in solution, and compute the contact maps of the two models, and plot the difference
of the maps:
```@example contacts
pdb = wget("2cpb", "model 1 2")
models = collect(eachmodel(pdb))
c1 = contact_map(models[1])
c2 = contact_map(models[2])
c_diff = c2 - c1
heatmap(c_diff)
```

A difference of distance maps can be similarly obtained by computing continous contact maps:
```@example contacts
c1 = contact_map(models[1]; discrete=false)
c2 = contact_map(models[2]; discrete=false)
c_diff = c2 - c1
heatmap(c_diff)
```

## Customizing the plot

All `heatmap` parametes can be customized using the `Plots` keyword syntax. Above, 
we illustrated this by removing the color bar and changing the color scale. 

Common customization options are:

- `xstep`: the stride of the x-axis ticks 
- `ystep`: the stride of the y-axis ticks
- `color`: the color palette to use (default: `:grayC` for distances, `:Greys_9` for binary maps)
- `clims`: the range of the color scale.
- `colorbar_title`: the title of the colorbar. Default: "distance (Å)" for distances, no title for binary maps.

## Indexing

The `ContactMap` data structure can be indexed to extract the contacts of a specific 
residue. For example:

```@example
using PDBTools
ats = read_pdb(PDBTools.DIMERPDB);
cA = select(ats, "chain A");
cB = select(ats, "chain B");
map = contact_map(cA, cB; discrete=false, dmax=12.0)
map[235,:] # all distances below 12.0 Angs of residue 235 of cA with cB
```

## Data structure and auxiliary functions

```@docs
ContactMap
residue_residue_distance
```
