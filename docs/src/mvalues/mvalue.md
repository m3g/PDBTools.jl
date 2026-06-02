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

Here we implement the Moeser/Horinek for urea, Auton/Bolen, and a version of Moeser/Horinek without Gly-non-ideality corrections for other cosolvents ([see this section](@ref mh_app)) models.
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
native_state = read_pdb(PDBTools.src_dir*"/tools/mvalue/testing/1MJC_native.pdb", "protein")
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
