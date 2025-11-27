# Structure object

```@docs
Structure
```

## Definition of structure object

A `Structure` object can be used to store an array of `Atom` objects along with, for example,
the unitcell information or other metadata: 

```@example structure
using PDBTools
ats = read_pdb(PDBTools.TESTPBC)
uc = read_unitcell(PDBTools.TESTPBC)
str = Structure(ats; unitcell=uc)
```

## Indexing and iteration

A `Structure` object behaves like a regular `Vector{Atom}` in most contexts, as it implements the `AbstractVector` interface. For example, the `str` object above can be indexed and iterated as a regular vector of atoms.

To fetch atom by indexing, one can do:

```@example structure
str[begin] # or str[1]
```

```@example structure
str[10]
```

```@example structure
str[end]
```

And all iterators that apply to raw `Atom` vectors also apply to the `Structure` object. For example, let us collect the chains of the structure:

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





