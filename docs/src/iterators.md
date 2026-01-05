```@meta
CollapsedDocStrings = true
```

# Iterators

PDBTools.jl provides lazy iterators over Residues, Chains, Segments, and Models of a structure file. The iterators behave similarly, and can be used bo computed properties of independent structural elements. The documentation bellow exemplifies in more detail the features associated to Residue and Chain interators, but the properties and valid for Segment and Model iterators similarly.

## Iterate over residues (or molecules)

The `eachresidue` iterator enables iteration over the residues of a structure. In PDB files, distinct molecules are often treated as separate residues, so this iterator can be used to iterate over the molecules within a structure. For example:

```@example iterators
using PDBTools
protein = read_pdb(PDBTools.SMALLPDB);
count(res -> resname(res) == "ALA", eachresidue(protein))
```

Here, we use `eachresidue` to count the number of residues named "ALA". This highlights the distinction between residue-level and atom-level operations.

### Collecting Residues into a Vector

Residues produced by `eachresidue` can be collected into a vector for further processing:

```@example iterators
residues = collect(eachresidue(protein))
```

and the atoms of a specific residue can be seen by indexing the residue:

```@example iterators
residues[1]
```

!!! note
    Iterators or collected vectors *do not* create copies of the original atom data. This means that any changes made to the residue vector will directly modify the corresponding data in the original atom vector.

### Iterating Over Atoms Within Residues

You can iterate over the atoms of one or more residues using nested loops. Here, we compute the total number of atoms of ALA residues: 

```@example iterators
let n_ala_cys = 0
    for residue in eachresidue(protein)
        if name(residue) in ("ALA", "CYS")
            for atom in residue
                n_ala_cys += 1
            end
        end
    end
    n_ala_cys
end
```

This method produces the same result as the more concise approach:

```@example iterators
sum(length(r) for r in eachresidue(protein) if name(r) in ("ALA", "CYS"))
```

Alternativelly, an image (not a copy!) of the atoms corresponding to a residue can be obtained with `get_atoms`:

```@example iterators
r1_atoms = get_atoms(residues[1])
```

### Reference documentation

```@docs
Residue
eachresidue
```

#### Residue property getters

The following functions provide convenient access to properties of `Residue` objects:

```@docs
name(::Residue)
resname(::Residue)
residue(::Residue)
resnum(::Residue)
chain(::Residue)
model(::Residue)
segname(::Residue)
mass(::Residue)
get_atoms(::Residue)
charge(::Residue)
residuename
```

## Iterate over chains

The `eachchain` iterator in PDBTools enables users to iterate over the chains in a PDB structure. A PDB file may contain multiple protein chains. This iterator simplifies operations involving individual chains.

```@example iterators
ats = read_pdb(PDBTools.CHAINSPDB);
chain.(eachchain(ats)) # Retrieve the names of all chains in the structure
```

```@example iterators
model.(eachchain(ats)) # Retrieve the model numbers associated with each chain
```

```@example iterators
chain_A1 = first(eachchain(ats)); # Access the first chain in the iterator
```

```@example iterators
resname.(eachresidue(chain_A1)) # Retrieve residue names for chain A in model 1
```

In the example above, the `chain.` command retrieves the names of all chains in the structure, while  `model.` command lists the model numbers for each chain. This PDB structure contains two models for chain A, where the third residue changes from leucine (LEU) in model 1 to valine (VAL) in model 2.

### Collect chains and indexing

As seen in the previous example, The `first` and `last` commands allow quick access to the first an last elements in the iterator. For more specific indexing, you can collect all chains into an array and then use numerical indices to access them.

```@example iterators
chains = collect(eachchain(ats))
```

```@example iterators
chain_B = chains[2]
```

### Modifying Atom Properties in a Chain

Any changes made to the atoms of a chain variable directly overwrite the properties of the original atoms in the structure. For example, modifying the occupancy and beta-factor columns of atoms in model 2 of chain A will update the corresponding properties in the original structure.

In the example below, the `occup` and `beta` properties of all atoms in model 2 of chain A are set to 0.00. The changes are reflected in the original `ats` vector, demonstrating that the modifications propagate to the parent data structure.

```@example iterators
first(eachchain(ats))
```

```@example iterators
for chain in eachchain(ats)
    if name(chain) == "A" && model(chain) == 1
        for atom in chain
            atom.occup = 0.00
            atom.beta = 0.00
        end
    end
end
first(eachchain(ats))
```

This behavior ensures efficient data manipulation but requires careful handling to avoid unintended changes. 

### Reference documentation

```@docs
Chain
eachchain
```

#### Chain property getters

The following functions provide convenient access to properties of `Chain` objects:

```@docs
name(::Chain)
chain(::Chain)
model(::Chain)
segname(::Chain)
mass(::Chain)
get_atoms(::Chain)
```

## Iterate over segments 

The `eachsegment` iterator enables iteration over the segments of a structure. For example:

```@example iterators
read_pdb(PDBTools.DIMERPDB)
eachsegment(ats)
```

```@example iterators
name.(eachsegment(ats))
```

The result of the iterator can also be collected, with:
```@example iterators
s = collect(eachsegment(ats))
```

```@example iterators
s[1]
```

These segment structure *does not* copy the data from the original atom vector. Therefore, changes performed on these vectors will be reflected on the original data.  

Iterators can be used to obtain or modify properties of the segments. Here we illustrate computing the mass of
each segment and renaming segment of all atoms with the segment indices:

```@example iterators
s = collect(eachsegment(ats))
```

Properties of each segment can then be obtained by broadcasting over the segments:

```@example iterators
mass.(s)
```

```@example iterators
formula.(s)
```

And iterating over the segments can allow changing properties of the atoms in a segment-specific way. For instance, here we change the segment names:

```@example iterators
for (iseg, seg) in enumerate(eachsegment(ats))
    for at in seg
        at.segname = "$(at.segname)$iseg"
    end
end
```

```@example iterators
collect(eachsegment(ats))
```

### Reference documentation

```@docs
Segment
eachsegment
```

#### Segment property getters

The following functions provide convenient access to properties of `Segment` objects:

```@docs
name(::Segment)
segname(::Segment)
mass(::Segment)
get_atoms(::Segment)
```

## Iterate over models

The `eachmodel` iterator enables iteration over the segments of a structure. For example:

```@example iterators
ats = wget("8S8N");
eachmodel(ats)
```

```@example iterators
model.(eachmodel(ats))
```

The result of the iterator can also be collected, with:

```@example iterators
m = collect(eachmodel(ats))
```

```@example iterators
m[1]
```

The model structure *does not* copy the data from the original atom vector. Therefore, changes performed on these vectors will be reflected on the original data.  

Iterators can be used to obtain or modify properties of the models. Here we illustrate computing the mass of
each segment and renaming segment of all atoms with the segment indices:

```@example iterators
center_of_mass.(eachmodel(ats))
```

### Reference documentation

```@docs
Model
eachmodel
```

#### Model property getters

The following functions provide convenient access to properties of `Model` objects:

```@docs
model(::Model)
get_atoms(::Model)
```

