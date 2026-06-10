```@meta
CollapsedDocStrings = true
```

# Protein conformational changes

Consider these two states of a model protein, a native and a denatured (straight chain) state, obtained from  a simulation. Here the conformational change could be of any kind. We load the structures of the two states:
```@example mvalue
using PDBTools
native_state = read_pdb(PDBTools.MJC_NATIVE, "protein")
desnat_state = read_pdb(PDBTools.MJC_DESNAT, "protein")
```

The denatured state has a greater surface area than the native state. Thus, cosolvents 
that bind preferentially to the surface, as urea, should promote a stabilization of the
denatured state. This is obtained with:
```@example mvalue
m = mvalue(native_state, desnat_state, "urea"; model=MoeserHorinek)
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
<img src="../../assets/mvalue.png" width=30%>
```
