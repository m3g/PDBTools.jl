# Structure object

```@docs
Structure
```

## Definition of structure object

A `Structure` object can be used to store an array of `Atom` objects along with, for example,
the unitcell information or other metadata. A `Structure` object behaves like a regular
`Vector{Atom}` in most contexts, as it implements the `AbstractVector` interface:

```@example structure
using PDBTools
ats = read_pdb(PDBTools.TESTPBC)
uc = read_unitcell(PDBTools.TESTPBC)
str = Structure(ats; unitcell=uc)
```

## Indexing and iteration

The `str` object above can be iterated and indexed as regular vector of atoms:

```@example structure
str[1]
```

```@example structure
str[end]
```

```@example structure
collect(eachchain(str))
```

## Data field assignment and access

Additionally, the `unitcell` field (or any other property defined), can be obtained by directly accessing the field
with the corresponding name:

```@example structure
str.unitcell
```

More fields can be added to the `Structure` object by simple assignment:

```@example structure
str.filename = "file.pdb"
```

which can instantly be acessed by:

```@example structure
str.filename
```





