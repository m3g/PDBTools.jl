```@meta
CollapsedDocStrings = true
```

# Tools

These tools may call external programs to perform each task. Please verify the installation of 
the necessary tool for each case. 

## Add hydrogens with OpenBabel

```@docs
add_hydrogens!
```

## Custom protein residue types

It is possible to add to the list of protein residues, custom residue types. 
This can be done by simply adding to the `PDBTools.protein_residues` dictionary
of residues a new `PDBTools.ProteinResidue` entry. 

```@docs
add_protein_residue!
remove_custom_protein_residues!
```

For example, here we create
a new resiude type `NEW` with the same properties of an `ALA` residue. To 
remove all custom protein residues, use `remove_custom_protein_residues!()`.

```@example custom_types
using PDBTools
add_protein_residue!("NEW", PDBTools.protein_residues["ALA"])
```

Then, an atom of residue name `NEW` will be recognized as a protein residue:
```@example custom_types
atom = Atom(resname="NEW")
isprotein(atom)
```

To remove the elements added, do:

```@example custom_types
remove_custom_protein_residues!()
```

## The SIRAH force-field residues and element types

Conveniencie functions can be created to add sets of new types of residues and atom types
to the list of residues and elements. This is illustrated in the 
[`custom_types.jl`](https://github.com/m3g/PDBTools.jl/blob/main/src/custom_types.jl) file of the source code, in this case for the residues and atom
types of the [`SIRAH`](http://www.sirahff.com/) force field for Coarse-Grained protein simulations.

!!! note
    - Masses of the SIRAH beads are not useful in general, and do not have a spcific physical meaning. 
    - Residue `sX` is interpreted as a bridged Cysteine.
    - To compute SASAs of SIRAH models, be careful in adding the SIRAH custom elements first, and use
      the `sasa_particles(SIRAH, ...)` method.

### Defining atom and residue types

Here we repeteadly call `remove_custom_residues!()` and `remove_custom_elements!()` to guarantee the proper execution of the
test codes, this is not necessary in a regular use of these functions.

With those definitions, adding all SIRAH protein residue types and element names can be done with:
```@example sirah
using PDBTools 
custom_protein_residues!(SIRAH)
custom_elements!(SIRAH)
```

With that, SIRAH structures have a special support. Let us start reading a CG model:

```@example sirah
sirah_pdb = read_pdb(PDBTools.SIRAHPDB)
resname.(eachresidue(sirah_pdb))
```

Note that the residue names of the SIRAH force-field are non-standard (`sI`, `sR`, etc.), but the sequence
is properly retrieved with standard one-letter codes: 

```@example sirah
getseq(sirah_pdb)
```

And all protein atoms are recognized:

```@example sirah
all(isprotein.(sirah_pdb))
```

### Computing SIRAH solvent accessible area

To compute the solvent accessible surface area of SIRAH models, call the `sasa_particles(SIRAH, ...)` method:

```@example sirah
s_sirah = sasa_particles(SIRAH, sirah_pdb)
```

which can be decomposed as usual:
```@example sirah
sasa(s_sirah, "sidechain")
```

### Removing the custom types

To remove the custom element types and residues, do:

```@example sirah
remove_custom_protein_residues!()
remove_custom_elements!()
```

## Move atoms and center of mass

The `center_of_mass` function can be used to compute the center of mass of set of atoms, and the 
`moveto!` function can be used to move the center of mass of the atoms to the origin (by default) 
or to a specified position:

```@docs
center_of_mass
moveto!
```
