```@meta
CollapsedDocStrings = true
```

# [Protein transfer free energy and m-values](@id mvalues)

These function compute the transfer free energies of proteins from water to different solvent, using 
the Tanford additive transfer models. 

*m*-values are the transfer free-energy difference, here in `kcal/mol`, between two structures from water
to a 1M solution of a cosolvent. Estimates of `m-values` for denaturation events can be computed
using the Creamer estimates for denatued accessible surface areas.

```@docs
transfer_free_energy
mvalue
```

Here we implement four additive transfer models: the established Auton/Bolen model, the universal-backbone
Moeser/Horinek model, and the
`Accessibility` model, which explicitly accounts for the mutual shielding between backbone and side-chain
groups ([see this section](@ref accessibility_model)).
([1](https://doi.org/10.1021/jp409934q), 
[2](https://doi.org/10.1016/s0076-6879(07)28023-1), 
[3](https://www.pnas.org/doi/10.1073/pnas.0706251104)). 
Typically, these models are used to obtain the effect of cosolvent on the structural stability of proteins,
but the current implementation allows the practical use of these functions to compute *m*-values of 
more general transformations, as described in the examples.

The transfer free energy of a protein from water to a 1M solution of a cosolvent can be estimated
with the `transfer_free_energy` function:

```@example mvalue
using PDBTools
native_state = read_pdb(PDBTools.MJC_NATIVE, "protein")
tfe = transfer_free_energy(native_state, "urea")
```

The resulting `TransferFreeEnegy` object contains the information of the contribution of each residue
to the transfer free energy obtained, split into backbone and side-chain contributions:
```@example mvalue
tfe.residue_contributions_bb[1]
```
```@example mvalue
tfe.residue_contributions_sc[1]
```

When multiple protein conformations are of interest, the `mvalue` methods provide a direct way to compute the 
variations in transfer free energies associated to the states involved, as shown in the following examples.

## Save and load TFE/m-value data

`TransferFreeEnergy` or `MValue` objects can be saved to a json file with `save` and restored with `load`:

```@example mvalue
outfile = tempname() * ".json"
save(outfile, tfe)
tfe_loaded = load(TransferFreeEnergy, outfile)
```

The saved file includes the transfer model name (for example, `AutonBolen` or `MoeserHorinek`),
which is used when loading to reconstruct the original parametric type:

```@example mvalue
typeof(tfe)
```

```@example mvalue
typeof(tfe_loaded)
```

The same interface can be used for `MValue` objects:

```@example mvalue
desnat_state = read_pdb(PDBTools.MJC_DESNAT, "protein")
m = mvalue(native_state, desnat_state, "urea")
save(outfile, m)
m_loaded = load(MValue, outfile)
```
