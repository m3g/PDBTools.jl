```@meta
CollapsedDocStrings = true
```

# [Secondary structure](@id secondary-structure)

These functions provide an interface to compute the secondary structure assignment of proteins
using the STRIDE and DSSP algorithms, using as inputs vectors of `PDBTools.Atom`s.

```@docs
stride_run
dssp_run
```

The `stride_run` and `dssp_run` functions return a vector of `SSData` objects, each containing the secondary structure assignment and backbone dihedral angles for each residue.

These functions return a vector of `SSData` objects, as defined in [ProteinSecondaryStructures.jl](https://biojulia.dev/ProteinSecondaryStructures.jl). The secondary structure assignment codes are available 
in the [ProteinSecondaryStructures.jl documentation](https://biojulia.dev/ProteinSecondaryStructures.jl/stable/overview/#Secondary-structure-classes).

!!! note
    Non-standard residue or atom names may lead to incorrect secondary structure assignments. To ensure accurate results, it is recommended to replace non-standard names with their standard three-letter codes before running these functions, and to verify that all backbone atoms (N, CA, C, O) are present in each residue. Errors or warnings will be issued if residue names are not recognized or backbone atoms are missing.

## Example using STRIDE

The `stride_run` function runs the STRIDE algorithm on the provided array of atoms.
```@example secondary-structure
using PDBTools
atoms = read_pdb(PDBTools.TESTPDB, "protein")
ss = stride_run(atoms)
```

To run with DSSP just use the `dssp_run` function instead. 

## Utility functions 

`PDBTools` also reexports the `ss_composition` , `ss_name`, `ss_code`, `ss_number` of `ProteinSecondaryStructures.jl`, which can be used to analyze the secondary structure assignment results. For example, to compute the secondary structure composition from the `ss` vector obtained above:

```@example secondary-structure
ss_composition(ss)
```

```@example secondary-structure
ss_name(ss[1]) # name of the secondary structure of the first residue
```

For further information refer to to [ProteinSecondaryStructures.jl documentation](https://biojulia.dev/ProteinSecondaryStructures.jl/stable/user_guide/#Retrieving-names,-codes,-and-numeric-codes).




