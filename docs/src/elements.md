```@meta
CollapsedDocStrings = true
```

# Atomic and molecular properties

## Atom properties

Some simple atom properties can be retrieved using special functions, which
operate on atoms of the type `Atom`. For example:

```julia-repl
julia> atoms = read_pdb("./file.pdb");

julia> printatom(atoms[1])
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1   OW     SOL     X        1        1   54.370   45.310   33.970  0.00  0.00     1       -         1

julia> mass(atoms[1])
14.0067

julia> atomic_number(atoms[1])
7

julia> element(atoms[1])
"N"

julia> element_name(atoms[1])
"Nitrogen"
```

```@docs
mass
element
element_name
element_symbol
element_symbol_string
element_vdw_radius
printatom
```

## Molecular formula and stoichiometry

The formula or stoichiometry of a selection can be obtained from a vector of Atoms:

```julia-repl
julia> atoms = wget("1LBD","protein and residue 1");

julia> f = formula(atoms)
C₃N₁O₂

julia> stoichiometry(select(atoms,"water"))
H₂O₁

```

```@docs
formula
stoichiometry
```

## [Custom Atom fields](@id custom-atom-fields)

Custom atom fields can be created in `Atom` objects by defining the `custom` keyword.
By default, `custom == nothing`. The custom fields can be added on construction, or 
with the `add_custom_field` function, which creates a new instance of an `Atom` 
with the added value in the custom field:

For example:

```jldoctest
julia> using PDBTools

julia> atom = Atom(custom="TEST");

julia> atom.custom
"TEST"

julia> atom = Atom(;name = "CA", resname="ALA"); # no custom field

julia> atom.resname
"ALA"

julia> new_atom = add_custom_field(atom, Dict(:charge => 2.0));

julia> new_atom.resname
"ALA"

julia> new_atom.custom[:charge]
2.0

```

```@docs
add_custom_field
```

# Elements for custom atom types

The types of atoms that `PDBTools` recognizes is defined in the `PDBTools.elements` dictionary. 
If new atom types are defined, it is possible to add these types to the dictionary, such that
other functions work for the new types. The function to be used is `add_element!`.

```@docs
add_element!
remove_custom_elements!
```

## Additional property-retrieving and set functions

The following functions are supported as part of the API, as a intending to interface
with `AtomsBase`. Nevertheless, currently these functions do not overload the exported
ones from `AtomsBase`, because that package is in a unstable state.

| Function   |  Example              |  Output |
|:-----------|:----------------------|:-------:|
|`atomic_number(::Atom)` | `atomic_number(Atom(name="NE2"))` |  `7` |
|`atomic_symbol(::Atom)` |  `atomic_symbol(Atom(name="NE2"))` |  `:N` |
|`atomic_mass(::Atom)`   |  `atomic_mass(Atom(name="NE2"))` |  `14.0067` |

```@docs
atomic_number
atomic_symbol
atomic_mass
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

## [The SIRAH force-field residues and element types](@id sirah)

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

Here we repeteadly call `remove_custom_residues!()` and `remove_custom_elements!()` to guarantee the proper execution of the test codes, this is not necessary in a regular use of these functions.

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

### Removing the custom types

To remove the custom element types and residues, do:

```@example sirah
remove_custom_protein_residues!()
remove_custom_elements!()
```




