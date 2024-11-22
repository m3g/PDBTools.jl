# [Examples](@id examples)

## Selecting the active site of a protein

The `1BSX` pdb file is a structure that contains a dimer of the thyroid hormone
receptor-beta bound to the ligand T3. Here we select all residues of `chain A`,
which is one of the monomers, that within 3.5$\AA$ of the ligand:

```jldoctest
julia> using PDBTools

julia> atoms = wget("1BSX", "chain A"; format="PDB");

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
6-element Vector{InlineStrings.String7}:
 "PHE"
 "ARG"
 "LEU"
 "ASN"
 "LEU"
 "HIS"
```

!!! note
   - The `Atom[]` creates an empty vector of `PDBTools.Atom` objects, and we
     append to this array the list of atoms of each residue. 
   - We opt here to download the file in the `"PDB"` format, because the chain
     identifier in the `mmCIF` deposited file does not include the ligand in chain A.



## Storing partial charges

Here we exemplify the use of a custom field to store partial charges for all atoms in a protein:

```jldoctest
julia> using PDBTools

julia> ats = wget("1BSX", "protein");

julia> charges = ones(length(ats));

julia> ats_with_charges = add_custom_field.(ats, charges); # charges in custom field

julia> ats_with_charges[1].custom
1.0

```
