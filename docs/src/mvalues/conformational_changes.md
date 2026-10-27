```@meta
CollapsedDocStrings = true
```

# Protein conformational changes

Consider these two states of a model protein, a native and a denatured (straight chain) state, obtained from  a simulation. Here the conformational change could be of any kind. We load the structures of the two states:
```@example mvalue
using PDBTools, LazyArtifacts
native_state = read_pdb(PDBTools.test_file("1MJC_native.pdb"), "protein")
desnat_state = read_pdb(PDBTools.test_file("1MJC_straight.pdb"), "protein")
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

We can set the `beta` fields (for example) of the atoms as the residue contributions, and observe 
the contributions drawn on top of the structure:
```@example mvalue
for (ir, r) in enumerate(eachresidue(native_state)) # iterate over residues
    # total contribution of residue ir
    c_residue = m.residue_contributions_sc[ir] + m.residue_contributions_bb[ir]
    for at in r # iterate over atoms in residue
        at.beta = c_residue
    end
end
phe_and_trp = select(native_state, "resname PHE TRP")
visualize(
    native_state => (color_by=beta.(native_state),),
    phe_and_trp => (style=:ballandstick, color_by=beta.(phe_and_trp))
)
```
Where we can see that Phe53 is a buried hydrophobic residue with the most important contribution to the 
denaturation m-value. 