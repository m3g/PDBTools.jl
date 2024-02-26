# [Examples](@id examples)

## Selecting the active site of a protein

The `1BSX` pdb file is a structure that contains a dimer of the thyroid hormone
receptor-beta bound to the ligand T3. Here we select all residues of `chain A`,
which is one of the monomers, that within 3.5$\AA$ of the ligand:

```jldoctest
julia> using PDBTools

julia> atoms = wget("1BSX", "chain A");

julia> protein = select(atoms, "protein");

julia> ligand = select(atoms, "resname T3");

julia> active_site_atoms = Atom[]
       for residue in eachresidue(protein)
           if distance(residue, ligand) < 3.5
               append!(active_site_atoms, atom for atom in residue)
           end
       end

julia> length(active_site_atoms)
56

julia> resname.(eachresidue(active_site_atoms))
6-element Vector{String}:
 "PHE"
 "ARG"
 "LEU"
 "ASN"
 "LEU"
 "HIS"
```

Note that `Atom[]` creates an empty vector of `PDBTools.Atom` objects, and we
append to this array the list of atoms of each residue.

