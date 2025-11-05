```@meta
CollapsedDocStrings = true
```

# [*m*-value (protein transfer free energy) calculator](@id mvalues)

*m*-values are the transfer free-energy difference, here in `kcal/mol`, between two structures from water
to a 1M solution of a cosolvent.

```@docs
mvalue
```

The *m*-values can be estimated by the Tanford transfer model, in which each amino acid contributes
to the transfer free energy according to the change in its solvent accessible surface area and 
experimental values of individual amino acid transfer energies, and a backbone contribution
(here we implement the Auton/Bolen and Moeser/Horinek models 
[[1](https://doi.org/10.1021/acs.jpcb.7b02138), 
[2](https://doi.org/10.1016/s0076-6879(07)28023-1), 
[3](https://www.pnas.org/doi/10.1073/pnas.0706251104)]). 
Typically, these models are used to obtain the effect of cosolvent on the structural stability of proteins,
but the current implementation allows the practical use of these functions to compute *m*-values of 
more general transformations, as described in the examples.

## Protein denaturation

Consider these two states of a model protein, a native and a denatured (straight chain) state, for
which we compute the SASA:
```@example mvalue
using PDBTools
native_state = read_pdb(PDBTools.src_dir*"/tools/mvalue/1MJC_native.pdb", "protein")
sasa_native = sasa_particles(native_state)
```

```@example mvalue
desnat_state = read_pdb(PDBTools.src_dir*"/tools/mvalue/1MJC_straight.pdb", "protein")
sasa_desnat = sasa_particles(desnat_state)
```

The denatured state has a greater surface area than the native state. Thus, cosolvents 
that bind preferentially to the surface, as urea, should promote a stabilization of the
denatured state. This is obtained with:
```@example mvalue
m = mvalue(sasa_native, sasa_desnat, "urea"; model=MoeserHorinek)
```
Where the `tot`, `bb` and `sc` fields contain, respectively, the total, backbone and side-chain contributions.
The `MValue` object contains, additionally, the contribution of the side chain and backbone of 
each amino acid residue type for the *m*-value, in the `residue_contributions_bb` and `residue_contributions_sc` fields.

We can set the `beta` fields (for example) of the atoms as the residue contributions:
```@example mvalue
for (ir, r) in enumerate(eachresidue(native_state)) # iterate over residues
    # total contribution of residue ir
    c_residue = m.residue_contributions_sc[ir] + m.residue_contributions_bb[ir]
    for at in r # iterate over atoms in residue
        at.beta = c_residue
    end
end
write_pdb("contrib.pdb", native_state)
```

And with that get an image (here produced with VMD) of the contributions of the residues
to the transfer free energies:

```@raw html
<img src="../assets/mvalue.png" width=30%>
```

## Cofactor binding

The `1BSX` protein-data-bank structure contains a nuclear hormone receptor bound to a cofactor:
```@example mvalue
bsx = wget("1BSX")
collect(eachchain(bsx))
```
Chains `A` and `B` belong to the receptor and the cofactor. Let us understand the effect of cosolvents
on the association of these two chains. 

First, we compute the SASA of chains `A` and `B`, thus including chain `A` in bound to the cofactor:
```@example mvalue
cAB = sasa_particles(select(bsx, "chain A B"))
```
and, then, we compute the SASA of chain `A` without the cofactor:
```@example mvalue
cA_free = sasa_particles(select(bsx, "chain A"))
```
Note that the SASA of chain A in the bound state is smaller than that of the free state, as expected:
```@example mvalue
sasa(cAB, "chain A")
```

We now compute the $m$-value of **chain A** only, but using the surface areas computed in the two different states:
```@example mvalue 
mvalue(cAB, cA_free, "urea"; sel="chain A", model=MoeserHorinek)
```
where the `tot` field is negative, indicating that exposing the cofactor binding surface is slightly favorable in urea.  

The same applies to the cofactor, chain B:
```@example mvalue
cB_free = sasa_particles(select(bsx, "chain B"))
mvalue(cAB, cB_free, "urea"; sel="chain B", model=MoeserHorinek)
```
and the exposed surface of the cofactor is also slightly stabilized in urea. Urea, thus tends to
destabilize the binding of the cofactor to the receptor.

By contrast, in a cosolvent that tends to promote protein aggregation, we have:
```@example mvalue 
mvalue(cAB, cA_free, "Sucrose"; sel="chain A")
```
and
```@example mvalue 
mvalue(cAB, cB_free, "Sucrose"; sel="chain B")
```
and thus Sucrose can stabilize cofactor binding. We remark that the values obtained here
are very small, and this is intended to be only an illustrative example.

## Alternative SASA calculations

The following functions can be used to compute *m*-values from the variation of the SASA per residue 
type, which allow the use of external tools to compute the SASA. This is used mostly for testing 
purposes. The functions allow the use of SASAs obtained directly from the [Auton & Bolen server](http://best.bio.jhu.edu/mvalue/), or 
from Gromacs SASA calculations.

```@docs
PDBTools.mvalue_delta_sasa
PDBTools.delta_sasa_per_restype
PDBTools.parse_mvalue_server_sasa
PDBTools.gmx_delta_sasa_per_restype
```
