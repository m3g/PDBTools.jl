```@meta
CollapsedDocStrings = true
```

# Protein denaturation

```@docs
CreamerDenaturedModel
mvalue(::CreamerDenaturedModel, ::AbstractString)
```

The estimate of the m-value associated to protein denaturation (the change in transfer free energy when a protein undergoes denaturation) can be (and is usually) computed using the estimates of exposed residue surface area obtained by [Creamer](https://doi.org/10.1021/bi962819o). 

Here, we compute the m-value estimates using the Creamer model, but first wrapping the atom array in the `CreamerDenaturedModel` type, which defines also the denaturation extent to be considered, by default, "mean":

```@example mvalue
using PDBTools
prot = read_pdb(PDBTools.TESTPDB, "protein")
model = CreamerDenaturedModel(prot)
```
A second argument of the constructor defines the extent to the "minimal", "mean", "maximal", as described in the docstring below. 

The m-value of denaturation in a cosolvent can be computed, then, with:

```@example mvalue
m = mvalue(model, "urea"; model=AutonBolen)
println("m-value of denaturation = ", m.tot)
```
where the output is a dictionary containing the total, backbone, side-chain, and residue-type specific contributions to the transfer free energies, in `kcal/mol`. 