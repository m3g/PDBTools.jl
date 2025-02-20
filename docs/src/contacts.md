```@meta
CollapsedDocStrings = true
```
# Contacts and contact maps

```@docs
contact_map
```

The `contact_map` function computes a contact map for a structure or a pair of structures.
These structures are typically proteins, but any structures defined by sequences of residues
can be provided as inputs. 

A typical usage consists in computing the contact map and plotting it:

```@example
using PDBTools
using Plots
ats = read_pdb(PDBTools.DIMERPDB);
cA = select(ats, "chain A");
cB = select(ats, "chain B");
map = contact_map(cA, cB) # contact map between chains A and B
heatmap(map)
```

In the example above we opted to plot a discrete contact map, with the default contact 
distance `dmax=4.0`. Now we change two parameters: `discrete=false` and `dmax=12.0`, to
compute a distance map up to a greater distance:

```@example
using PDBTools
using Plots
ats = read_pdb(PDBTools.DIMERPDB);
cA = select(ats, "chain A");
cB = select(ats, "chain B");
map = contact_map(cA, cB; discrete=false, dmax=12.0) # contact map between chains A and B
heatmap(map)
```

Similarly, we can produce plots for the contact map of a single structure. Here, we 
showcase the use of the `gap` parameter, to ignore residues closer in the sequence
by less than 4 residues:

```@example
using PDBTools
using Plots
ats = read_pdb(PDBTools.DIMERPDB);
cA = select(ats, "chain A");
distance_map = contact_map(cA; gap=4, discrete=false, dmax=12.0) # chain A only
discrete_map = contact_map(cA; gap=4, discrete=true, dmax=12.0) # chain A only
plot(
    heatmap(distance_map; colorbar=nothing, color=:davos), 
    heatmap(discrete_map); 
    layout=(1,2), size=(800,500)
)
```

All `heatmap` parametes can be customized using the `Plots` keyword syntax. Above, 
we illustrate this by removing the color bar and changing the color scale. 

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
