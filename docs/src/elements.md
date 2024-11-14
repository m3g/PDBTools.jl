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

## AtomsBase compatibility

!!! compat
    This interface requires at least PDBTools version 0.14.2.

The following functions are supported as part of the API, to conform the `AtomsBase` interface:

| Function   |  Example              |  Output |
|:-----------|:----------------------|:-------:|
|`atomic_number(::PDBTools.Atom)` | `atomic_number(Atom(name="NE2"))` |  `7` |
|`atomic_symbol(::PDBTools.Atom)` |  `atomic_symbol(Atom(name="NE2"))` |  `:N` |
|`atomic_mass(::PDBTools.Atom)`   |  `atomic_mass(Atom(name="NE2"))` |  `14.0067` |
|`position(::PDBTools.Atom)`      |  `position(Atom(name="NE2"))` |  `SVector{3,Float64}(0,0,0)` |


## Custom Atom fields

!!! compat
    Custom field support was introduced on PDBTools version 0.14.3.

Custom atom fields can be added to an `Atom` object by defining the `custom` dictionary.
The fields can be accessed by the standard dot syntax if the field name does not clash 
with an existing `Atom` field, or by the `custom_field` getter function. 

For example:

```julia-repl
julia> atom = Atom(index = 0; custom=Dict(:c => "c", :index => 1))
       0    X     XXX     X        0        0    0.000    0.000    0.000  0.00  0.00     0    XXXX         0

julia> atom.c
"c"

julia> atom.index
0

julia> custom_field(atom, :index)
1
```

Setting new custom fields follow the standard Julia dictionary syntax:

```julia-repl
julia> atom.custom[:new] = "NEW"
"NEW"

julia> atom.new
"NEW"

julia> custom_field(atom, :new)
"NEW"
```

!!! compat 
    The following feature was introduced in PDBTools version 0.14.4.

If a custom field with the `:mass` key is added to the atom, the `mass` function returns the mass
set at that field: 

```jldoctest
julia> using PDBTools

julia> atom = Atom();

julia> atom.custom[:mass] = 10.0
10.0

julia> mass(atom)
10.0
```

# Elements for custom atom types

!!! compat
    The `add_element!` function was introduced in version 1.4.0.

The types of atoms that `PDBTools` recognizes is defined in the `PDBTools.elements` dictionary. 
If new atom types are defined, it is possible to add these types to the dictionary, such that
other functions work for the new types. The function to be used is `add_element!`.

```@docs
add_element!
remove_custom_elements!
```




