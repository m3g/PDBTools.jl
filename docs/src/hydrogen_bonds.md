```@meta
CollapsedDocStrings = true
```

# Hydrogen-bond analysis

This tool provides a fast and practical way to compute hydrogen bonds in structures (if the structures contain hydrogen atoms).

```@docs
hydrogen_bonds  
```

## The `HBonds` structure

The `hydrogen_bonds` function returns an `OrderedDict{String, HBonds}`, where each key corresponds
to the selection pair used (e.g., `"protein => protein"` or `"protein => resname SOL"`).

Each `HBonds` object is an iterable container where each element is a named tuple with fields:
- `D`: Index of the donor atom
- `H`: Index of the hydrogen atom
- `A`: Index of the acceptor atom
- `r`: Distance between donor and acceptor (Å)
- `ang`: Angle between H-D and A-D vectors (degrees)

## Basic Usage

```@example hbonds
using PDBTools
pdb = read_pdb(PDBTools.test_dir*"/hbonds.pdb", "model 1")
uc = read_unitcell(PDBTools.test_dir*"/hbonds.pdb")
hbs = hydrogen_bonds(pdb, 
    "protein", # alias for "protein" => "protein"
    "protein" => "water",
    unitcell=uc # optional
)
```

Access the hydrogen bonds for a specific selection:

```@example hbonds
hb_prot = hbs["protein => protein"]
```

## Iterating over hydrogen bonds

The `HBonds` structure supports standard Julia iteration. You can iterate directly over
all hydrogen bonds:

```@example hbonds
for bond in hb_prot
    println("Donor: $(bond.D), Hydrogen: $(bond.H), Acceptor: $(bond.A), Distance: $(round(bond.r, digits=2)) Å")
    break  # just show first one for brevity
end
```

Use `enumerate` to get both index and bond:

```@example hbonds
for (i, bond) in enumerate(hb_prot)
    if i <= 3
        println("Bond $i: D=$(bond.D) -> A=$(bond.A), r=$(round(bond.r, digits=2)) Å, angle=$(round(bond.ang, digits=1))°")
    end
end
```

Collect all bonds into a vector:

```@example hbonds
hbonds_prot = collect(hb_prot)
``` 

## Finding specific bonds with `findfirst` and `findall`

Use `findfirst` to locate the first hydrogen bond matching a condition:

```@example hbonds
idx = findfirst(bond -> bond.r < 2.6, hb_prot)
```

Find first H-bond involving a specific donor atom index

```@example hbonds
acceptor_index = 55 
idx = findfirst(bond -> bond.A == acceptor_index, hb_prot)
```

The corresponding hydrogen bond being:

```@example hbonds
hb_prot[6]
``` 

Find all hydrogen bonds with angle less than 15 degrees

Use `findall` to get indices of all matching hydrogen bonds:

```@example hbonds
indices = findall(bond -> bond.ang < 15, hb_prot)
println("Found $(length(indices)) H-bonds with angle < 15°")
```

Find all hydrogen bonds shorter than 2.8 Å

```@example hbonds
short_hbonds = findall(bond -> bond.r < 2.8, hb_prot)
```

These indices can be used to extract the list from a `HBonds` data structure:

```@example hbonds
hb_prot[short_hbonds]
```

and will be similar to a `filter` operation, as shown below. 

## Filtering

Use `filter` to extract bonds matching criteria:

```@example hbonds
# Get all strong hydrogen bonds (short distance, small angle)
strong_bonds = filter(bond -> bond.r < 2.8 && bond.ang < 20, hb_prot)
```


### Use comprehensions for flexible queries:

Extract distances and angles of all hydrogen bonds:

```@example hbonds
distances = [(bond.r, bond.ang) for bond in hb_prot]
```

## Multiple selections

Compute hydrogen bonds between different groups of atoms:

```@example hbonds
hbs = hydrogen_bonds(pdb, "protein", "protein" => "resname HOH"; unitcell=uc)
```

### Intra-protein hydrogen bonds

```@example hbonds
hbs["protein => protein"]
```

### Protein-water hydrogen bonds  
```@example hbonds
hbs["protein => resname HOH"]
```

!!! note
    Different selections in a pair (like `"protein" => "resname SOL"`) must not have overlapping atoms.
    An error will be thrown if the same atom appears in both selections.
