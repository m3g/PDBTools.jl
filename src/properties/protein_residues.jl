#
# Data for natural protein residues
#
@kwdef struct ProteinResidue
    name::String
    three_letter_code::String
    one_letter_code::String
    type::String
    polar::Bool
    hydrophobic::Bool
    mono_isotopic_mass::Float64
    mass::Float64
    charge::Int
    custom::Bool = true
end

#! format: off
const protein_residues = OrderedDict{String,ProteinResidue}(
    "ALA" => ProteinResidue("Alanine",       "ALA", "A", "Aliphatic",  false, true,   71.037114,  71.0779,  0, false),
    "ARG" => ProteinResidue("Arginine",      "ARG", "R", "Basic",      true,  false, 156.101111, 156.1857,  1, false),
    "ASN" => ProteinResidue("Asparagine",    "ASN", "N", "Amide",      true,  false, 114.042927, 114.1026,  0, false),
    "ASP" => ProteinResidue("Aspartic acid", "ASP", "D", "Acidic",     true,  false, 115.026943, 115.0874, -1, false),
    "CYS" => ProteinResidue("Cysteine",      "CYS", "C", "Sulfuric",   true,  false, 103.009185, 103.1429,  0, false),
    "GLN" => ProteinResidue("Glutamine",     "GLN", "Q", "Amide",      true,  false, 128.058578, 128.1292,  0, false),
    "GLU" => ProteinResidue("Glutamic acid", "GLU", "E", "Acidic",     true,  false, 129.042593, 129.1140, -1, false),
    "GLY" => ProteinResidue("Glycine",       "GLY", "G", "Aliphatic",  true,  false,  57.021464,  57.0513,  0, false),
    "HIS" => ProteinResidue("Histidine",     "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "ILE" => ProteinResidue("Isoleucine",    "ILE", "I", "Aliphatic",  false, true,  113.084064, 113.1576,  0, false),
    "LEU" => ProteinResidue("Leucine",       "LEU", "L", "Aliphatic",  false, true,  113.084064, 113.1576,  0, false),
    "LYS" => ProteinResidue("Lysine",        "LYS", "K", "Basic",      true,  false, 128.094963, 128.1723,  1, false),
    "MET" => ProteinResidue("Methionine",    "MET", "M", "Sulfuric",   false, true,  131.040485, 131.1961,  0, false),
    "PHE" => ProteinResidue("Phenylalanine", "PHE", "F", "Aromatic",   false, true,  147.068414, 147.1739,  0, false),
    "PRO" => ProteinResidue("Proline",       "PRO", "P", "Cyclic",     false, false,  97.052764,  97.1152,  0, false),
    "SER" => ProteinResidue("Serine",        "SER", "S", "Hydroxylic", true,  false,  87.032028, 87.07730,  0, false),
    "THR" => ProteinResidue("Threonine",     "THR", "T", "Hydroxylic", true,  false, 101.047679, 101.1039,  0, false),
    "TRP" => ProteinResidue("Tryptophan",    "TRP", "W", "Aromatic",   false, true,  186.079313, 186.2099,  0, false),
    "TYR" => ProteinResidue("Tyrosine",      "TYR", "Y", "Aromatic",   true,  false, 163.063320, 163.1733,  0, false),
    "VAL" => ProteinResidue("Valine",        "VAL", "V", "Aliphatic",  false, true,   99.068414,  99.1311,  0, false),
    # Alternate protonation states for CHARMM and AMBER
    "ASPP" => ProteinResidue("Aspartic acid (protonated)", "ASP", "D", "Acidic", true,  false, 115.026943, 115.0874, 0, false),
    "GLUP" => ProteinResidue("Glutamic acid (protonated)", "GLU", "E", "Acidic", true,  false, 129.042593, 129.1140, 0, false),
    "HSD"  => ProteinResidue("Histidine (D)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HSE"  => ProteinResidue("Histidine (E)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HSP"  => ProteinResidue("Histidine (doubly protonated)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  1, false), 
    "HID"  => ProteinResidue("Histidine (D)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HIE"  => ProteinResidue("Histidine (E)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HIP"  => ProteinResidue("Histidine (doubly protonated)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  1, false), 
)
#! format: on



"""
    add_protein_residue!(resname::AbstractString, reference_residue::PDBTools.ProteinResidue)

Function to add a custom protein residue to the list of protein residues. The function will
return the `ProteinResidue` object that was added. To remove all custom protein residues
use `remove_custom_protein_residues!()`.

# Example

```jldoctest
julia> using PDBTools

julia> remove_custom_protein_residues!();

julia> add_protein_residue!("sA", PDBTools.protein_residues["ALA"])
PDBTools.ProteinResidue("sA", "ALA", "A", "Aliphatic", false, true, 71.037114, 71.0779, 0, true)

julia> isprotein(Atom(resname="sA"))
true

julia> remove_custom_protein_residues!(); # clean up
```

Here we repeatedly call `remove_custom_residues!()` to guarantee the proper execution of the
test codes, without any custom residues in the list of protein residues.

"""
function add_protein_residue!(resname::AbstractString, reference_residue::PDBTools.ProteinResidue)
    if haskey(PDBTools.protein_residues, resname)
        @warn """\n
            Residue $resname already exists in the list of protein residues. Overwriting.

        """ _file = nothing _line = nothing
    end
    PDBTools.protein_residues[resname] = PDBTools.ProteinResidue(
        resname,
        reference_residue.three_letter_code,
        reference_residue.one_letter_code,
        reference_residue.type,
        reference_residue.polar,
        reference_residue.hydrophobic,
        reference_residue.mono_isotopic_mass,
        reference_residue.mass,
        reference_residue.charge,
        true
    )
    return PDBTools.protein_residues[resname]
end

"""
    remove_custom_protein_residues!()

Function to remove all custom protein residues from the list of protein residues.

# Example

```jldoctest
julia> using PDBTools

julia> remove_custom_protein_residues!(); # clean up

julia> add_protein_residue!("sA", PDBTools.protein_residues["ALA"])
PDBTools.ProteinResidue("sA", "ALA", "A", "Aliphatic", false, true, 71.037114, 71.0779, 0, true)

julia> isprotein(Atom(resname="sA"))
true

julia> remove_custom_protein_residues!();

julia> isprotein(Atom(resname="sA"))
false
```

Here we repeatedly call `remove_custom_residues!()` to guarantee the proper execution of the
test codes, without any custom residues in the list of protein residues.

"""
remove_custom_protein_residues!() = filter!(r -> !last(r).custom, protein_residues)

