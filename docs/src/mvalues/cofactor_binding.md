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
