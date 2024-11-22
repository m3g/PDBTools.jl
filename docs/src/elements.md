```@meta
CollapsedDocStrings = true
```

# Atomic and molecular properties

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

The formula or stoichiometry of a selection can also be retrieved:

```julia-repl
julia> atoms = wget("1LBD","protein and residue 1");

julia> f = formula(atoms)
C₃N₁O₂

julia> stoichiometry(select(atoms,"water"))
H₂O₁

```

```@docs
mass
element
element_name
element_symbol
element_symbol_string
formula
stoichiometry
printatom
```

## Custom Atom fields

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

## Additional property-retrieving functions

The following functions are supported as part of the API, as a intending to interface
with `AtomsBase`. Nevertheless, currently these functions do not overload the exported
ones from `AtomsBase`, because that package is in a unstable state.

| Function   |  Example              |  Output |
|:-----------|:----------------------|:-------:|
|`atomic_number(::PDBTools.Atom)` | `atomic_number(Atom(name="NE2"))` |  `7` |
|`atomic_symbol(::PDBTools.Atom)` |  `atomic_symbol(Atom(name="NE2"))` |  `:N` |
|`atomic_mass(::PDBTools.Atom)`   |  `atomic_mass(Atom(name="NE2"))` |  `14.0067` |
|`position(::PDBTools.Atom)`      |  `position(Atom(name="NE2"))` |  `SVector{3,Float64}(0,0,0)` |


```@docs
atomic_number
atomic_symbol
atomic_mass
position
```




