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

## Storing partial charges

Here we exemplify the use of a custom field to store partial charges for all atoms
in a protein:

```julia-repl
julia> using PDBTools

julia> pdb = wget("1BSX", "protein");

julia> charges = rand(length(pdb));

julia> for (i, atom) in enumerate(pdb)
           atom.custom[:charge] = charges[i]
       end

julia> pdb[1].custom[:charge]
0.09441681249467149

julia> custom_field(pdb[1], :charge) # alternative getter function
0.09441681249467149

julia> custom_field.(pdb, :charge) # broadcast to get all charges (with the dot syntax)
3994-element Vector{Float64}:
 0.09441681249467149
 0.17811534472805368
 ⋮
 0.8254040639975442
 0.6153943592336552
```
