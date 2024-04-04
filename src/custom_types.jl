#
# Custom protein residues that should not be loaded by default
#

#
# Sirah force field protein residues and elements
#
struct SIRAH end

function custom_protein_residues!(::Type{SIRAH}; protein_residues=PDBTools.protein_residues)
    @warn """\n
        Residue `sX` will be interpreted as bridged Cysteine.

    """ _file = nothing _line = nothing
#! format: off
    protein_residues["sA"] = ProteinResidue("Alanine",       "ALA", "A", "Aliphatic",  false, false,  71.037114,  71.0779,  0)
    protein_residues["sR"] = ProteinResidue("Arginine",      "ARG", "R", "Basic",      true,  false, 156.101111, 156.1857,  1)
    protein_residues["sN"] = ProteinResidue("Asparagine",    "ASN", "N", "Amide",      true,  false, 114.042927, 114.1026,  0)
    protein_residues["sD"] = ProteinResidue("Aspartic acid", "ASP", "D", "Acidic",     true,  false, 115.026943, 115.0874, -1)
    protein_residues["sC"] = ProteinResidue("Cysteine",      "CYS", "C", "Sulfuric",   false, false, 103.009185, 103.1429,  0)
    protein_residues["sQ"] = ProteinResidue("Glutamine",     "GLN", "Q", "Amide",      true,  false, 128.058578, 128.1292,  0)
    protein_residues["sE"] = ProteinResidue("Glutamic acid", "GLU", "E", "Acidic",     true,  false, 129.042593, 129.1140, -1)
    protein_residues["sG"] = ProteinResidue("Glycine",       "GLY", "G", "Aliphatic",  false, false,  57.021464,  57.0513,  0)
    protein_residues["sHE"] = ProteinResidue("Histidine",     "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0)
    protein_residues["sI"] = ProteinResidue("Isoleucine",    "ILE", "I", "Aliphatic",  false, true,  113.084064, 113.1576,  0)
    protein_residues["sL"] = ProteinResidue("Leucine",       "LEU", "L", "Aliphatic",  false, true,  113.084064, 113.1576,  0)
    protein_residues["sK"] = ProteinResidue("Lysine",        "LYS", "K", "Basic",      true,  false, 128.094963, 128.1723,  1)
    protein_residues["sM"] = ProteinResidue("Methionine",    "MET", "M", "Sulfuric",   false, false, 131.040485, 131.1961,  0)
    protein_residues["sF"] = ProteinResidue("Phenylalanine", "PHE", "F", "Aromatic",   false, true,  147.068414, 147.1739,  0)
    protein_residues["sP"] = ProteinResidue("Proline",       "PRO", "P", "Cyclic",     false, false,  97.052764,  97.1152,  0)
    protein_residues["sS"] = ProteinResidue("Serine",        "SER", "S", "Hydroxylic", true,  false,  87.032028, 87.07730,  0)
    protein_residues["sT"] = ProteinResidue("Threonine",     "THR", "T", "Hydroxylic", true,  false, 101.047679, 101.1039,  0)
    protein_residues["sW"] = ProteinResidue("Tryptophan",    "TRP", "W", "Aromatic",   false, true,  186.079313, 186.2099,  0)
    protein_residues["sY"] = ProteinResidue("Tyrosine",      "TYR", "Y", "Aromatic",   true,  false, 163.063320, 163.1733,  0)
    protein_residues["sV"] = ProteinResidue("Valine",        "VAL", "V", "Aliphatic",  false, true,   99.068414,  99.1311,  0)
    protein_residues["sX"] = ProteinResidue("Cysteine - bridged",     "CYS", "C", "Sulfuric",   false, false, 103.009185, 103.1429,  0)
#! format: on
    return nothing
end

function custom_elements!(::Type{SIRAH}; elements=PDBTools.elements)
    add_element!("GN", elements["N"])
    add_element!("GC", elements["C"])
    add_element!("BCG", elements["C"])
    add_element!("GO", elements["O"])
    add_element!("BCZ", elements["C"])
    add_element!("BNN1", elements["N"])
    add_element!("BNN2", elements["N"])
    add_element!("BSG", elements["S"])
    add_element!("BCB", elements["C"])
    add_element!("BOG", elements["O"])
    add_element!("BPG", elements["H"])
    add_element!("BND", elements["N"])
    add_element!("BOD", elements["O"])
    add_element!("BOE1", elements["O"])
    add_element!("BOE2", elements["O"])
    add_element!("BCE1", elements["C"])
    add_element!("BCE2", elements["C"])
    add_element!("BCD", elements["C"])
    add_element!("BSD", elements["S"])
    add_element!("BNE", elements["N"])
    add_element!("BPE", elements["H"])
    add_element!("BCE", elements["C"])
    nothing
end

@testitem "SIRAH" begin
    using PDBTools
    pdb_file = PDBTools.SIRAHPDB
    custom_protein_residues!(SIRAH)
    custom_elements!(SIRAH)
    pdb = readPDB(pdb_file)
    @test length(isprotein.(pdb)) == 22
    @test element(pdb[1]) == "N"
    @test element(pdb[2]) == "C"
    @test element(pdb[3]) == "C"
end
