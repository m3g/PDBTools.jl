```@meta
CollapsedDocStrings = true
```

# Cofactor binding

The `1BSX` protein-data-bank structure contains a nuclear hormone receptor bound to a cofactor:
```@example mvalue
using PDBTools
bsx = wget("1BSX")
collect(eachchain(bsx))
```
Chains `A` and `B` belong to the receptor and the cofactor. Let us understand the effect of cosolvents
on the association of these two chains. 

First, we select chains `A` and `B`, thus including chain `A` bound to the cofactor, and computes 
its transfer free energy to urea:
```@example mvalue
cAB = select(bsx, "chain A B")
cAB_tfe = transfer_free_energy(cAB, "urea")
```

and, then, we select chain `A` without the cofactor:
```@example mvalue
cA_free = select(bsx, "chain A")
cA_tfe = transfer_free_energy(cA_free, "urea")
```
and the cofactor,
```@example mvalue
cB_free = select(bsx, "chain B")
cB_tfe = transfer_free_energy(cB_free, "urea")
```

The free energy associated with dissociation is then:
```@example mvalue
(cA_tfe.tot + cB_tfe.tot) - cAB_tfe.tot
```
The small negative value indicates that dissociation is favored by urea.

By contrast, in a cosolvent that tends to promote protein aggregation, we have:
```@example mvalue 
cAB_tfe = transfer_free_energy(cAB, "sucrose")
cA_tfe = transfer_free_energy(cA_free, "sucrose")
cB_tfe = transfer_free_energy(cB_free, "sucrose")
(cA_tfe.tot + cB_tfe.tot) - cAB_tfe.tot
```
and thus Sucrose can stabilize cofactor binding. We remark that the values obtained here
are very small, and this is intended to be only an illustrative example.
