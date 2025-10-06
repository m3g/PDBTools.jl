```@meta
CollapsedDocStrings = true
```

# Tools

These tools may call external programs to perform each task. Please verify the installation of 
the necessary tool for each case. 

## Solvent Accessible Surface Area (SASA)

```@docs
sasa_particles
sasa
```

A typical run of these functions consists in providing the structure of a protein to the first function, `sasa_particles`, to obtain a `SASA` object, which contains the accessible area per atom:

```@example sasa
using PDBTools
prot = read_pdb(PDBTools.TESTPDB, "protein")
atom_sasa = sasa_particles(prot)
```

The `atom_sasa` object created above can be used to extract the total accessible area or the accessible area of any sub-surface. The `sasa` function provides an interface for those extractions:

```@example sasa
sasa(atom_sasa) # total
```
```@example sasa
sasa(atom_sasa, "polar") 
```

```@example sasa
sasa(atom_sasa, "backbone")
```

```@example sasa
sasa(atom_sasa, "resname THR and residue < 50") 
```

In some situations, it might be useful to visualize the surface. The dots that form the surface can be obtained by running `sasa_particles` with the `output_dots` option set to `true`. Here, we use fewer dots for better visualization:

```@example sasa
atom_sasa = sasa_particles(prot; n_dots=100, output_dots=true) 
```

Where the `atom_sasa.dots` field contais the dots that are accessible to the surface for each atom. These can be plotted, for example, with:
```@example sasa
using Plots
dots = reduce(vcat, atom_sasa.dots)
scatter(Tuple.(coor.(prot)); color=:orange, msw=0, label="") # atom coordinates
scatter!(Tuple.(dots); # surface dots
    color=:blue, ms=1, msw=0, ma=0.5, # marker properties
    label="",
)
```

## Check the stereochemistry of protein residues

```@docs
zeta
zeta_check
```

## Add hydrogens with OpenBabel

```@docs
add_hydrogens!
```

## Custom protein residue types

It is possible to add to the list of protein residues, custom residue types. 
This can be done by simply adding to the `PDBTools.protein_residues` dictionary
of residues a new `PDBTools.ProteinResidue` entry. For example, here we create
a new resiude type `NEW` with the same properties of an `ALA` residue. To 
remove all custom protein residues, use `remove_custom_protein_residues!()`.

```jldoctest
julia> using PDBTools

julia> remove_custom_protein_residues!();

julia> add_protein_residue!("NEW", PDBTools.protein_residues["ALA"])
PDBTools.ProteinResidue("NEW", "ALA", "A", "Aliphatic", false, true, 71.037114, 71.0779, 0, true)

julia> atom = Atom(resname="NEW");

julia> isprotein(atom)
true

julia> remove_custom_protein_residues!();
```

Here we repeteadly call `remove_custom_residues!()` to guarantee the proper execution of the
test codes, without any custom residues in the list of protein residues.

```@docs
add_protein_residue!
remove_custom_protein_residues!
```

### The SIRAH force-field residues and element types

Conveniencie functions can be created to add sets of new types of residues and atom types
to the list of residues and elements. This is illustrated in the 
[`custom_types.jl`](https://github.com/m3g/PDBTools.jl/blob/main/src/custom_types.jl) file of the source code, in this case for the residues and atom
types of the [`SIRAH`](http://www.sirahff.com/) force field for Coarse-Grained protein simulations.

With those definitions, adding all SIRAH protein residue types and element names can be done with:
```jldoctest
julia> using PDBTools 

julia> remove_custom_protein_residues!(); remove_custom_elements!();

julia> custom_protein_residues!(SIRAH)
┌ Warning: 
│ 
│     Residue `sX` will be interpreted as bridged Cysteine.
│ 
└ @ PDBTools

julia> custom_elements!(SIRAH)
┌ Warning:
│
│     The element masses are not the coarse-grained ones. This must be fixed in the future.
│
└ @ PDBTools

julia> sirah_pdb = read_pdb(PDBTools.SIRAHPDB);

julia> resname.(eachresidue(sirah_pdb))
5-element Vector{InlineStrings.String7}:
 "sI"
 "sR"
 "sX"
 "sI"
 "sG"

julia> getseq(sirah_pdb)
5-element Vector{String}:
 "I"
 "R"
 "C"
 "I"
 "G"

julia> all(isprotein.(sirah_pdb))
true

julia> remove_custom_protein_residues!(); remove_custom_elements!();
```

Note that the residue names of the SIRAH force-field are non-standard (`sI`, `sR`, etc.), but the sequence
is properly retrieved with standard one-letter codes, and all the atoms of the structure are recognized 
as being "protein" atoms.

Here we repeteadly call `remove_custom_residues!()` and `remove_custom_elements!()` to guarantee the proper execution of the
test codes.

## Move atoms and center of mass

The `center_of_mass` function can be used to compute the center of mass of set of atoms, and the 
`moveto!` function can be used to move the center of mass of the atoms to the origin (by default) 
or to a specified position:

```@docs
center_of_mass
moveto!
```
