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
    add_protein_residue!("sA",PDBTools.protein_residues["ALA"])
    add_protein_residue!("sR",PDBTools.protein_residues["ARG"])
    add_protein_residue!("sN",PDBTools.protein_residues["ASN"])
    add_protein_residue!("sD",PDBTools.protein_residues["ASP"])
    add_protein_residue!("sC",PDBTools.protein_residues["CYS"])
    add_protein_residue!("sQ",PDBTools.protein_residues["GLN"])
    add_protein_residue!("sE",PDBTools.protein_residues["GLU"])
    add_protein_residue!("sG",PDBTools.protein_residues["GLY"])
    add_protein_residue!("sHe",PDBTools.protein_residues["HIS"])
    add_protein_residue!("sI",PDBTools.protein_residues["ILE"])
    add_protein_residue!("sL",PDBTools.protein_residues["LEU"])
    add_protein_residue!("sK",PDBTools.protein_residues["LYS"])
    add_protein_residue!("sM",PDBTools.protein_residues["MET"])
    add_protein_residue!("sF",PDBTools.protein_residues["PHE"])
    add_protein_residue!("sP",PDBTools.protein_residues["PRO"])
    add_protein_residue!("sS",PDBTools.protein_residues["SER"])
    add_protein_residue!("sT",PDBTools.protein_residues["THR"])
    add_protein_residue!("sW",PDBTools.protein_residues["TRP"])
    add_protein_residue!("sY",PDBTools.protein_residues["TYR"])
    add_protein_residue!("sV",PDBTools.protein_residues["VAL"])
    add_protein_residue!("sX",
        PDBTools.ProteinResidue(
            name = "Cysteine - bridged",     
            three_letter_code = "CYS", 
            one_letter_code = "C", 
            type = "Sulfuric",   
            polar = false, 
            hydrophobic = false, 
            mono_isotopic_mass = 103.009185, 
            mass = 103.1429,  
            charge = 0, 
            custom = true
        )
    )
#! format: on
    return nothing
end

function custom_elements!(::Type{SIRAH}; elements=PDBTools.elements)
    @warn """\n
        The element masses are not the coarse-grained ones. This must be fixed in the future.

    """ _file = nothing _line = nothing
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
    add_element!("LN", elements["N"])
    add_element!("LN1", elements["N"])
    add_element!("BC12", elements["C"])
    add_element!("BC13", elements["C"])
    add_element!("BC14", elements["C"])
    add_element!("BC15", elements["C"])
    add_element!("BCT1", elements["C"])
    add_element!("BC21", elements["C"])
    add_element!("BC22", elements["C"])
    add_element!("BC23", elements["C"])
    add_element!("BC24", elements["C"])
    add_element!("BCT2", elements["C"])
    add_element!("BFO1", elements["O"])
    add_element!("BFO2", elements["O"])
    add_element!("WN1", elements["N"])
    add_element!("WN2", elements["N"])
    add_element!("WP1", elements["H"])
    add_element!("WP2", elements["H"])
    add_element!("NaW", elements["Na"])
    add_element!("ClW", elements["Cl"])
    add_element!("LN2", elements["N"])
    add_element!("LP1", elements["H"])
    add_element!("LP2", elements["H"])
    nothing
end

@testitem "SIRAH" begin
    using PDBTools
    pdb_file = PDBTools.SIRAHPDB
    remove_custom_protein_residues!()
    remove_custom_elements!()
    custom_protein_residues!(SIRAH)
    custom_elements!(SIRAH)
    pdb = readPDB(pdb_file)
    @test length(isprotein.(pdb)) == 22
    @test element(pdb[1]) == "N"
    @test element(pdb[2]) == "C"
    @test element(pdb[3]) == "C"
    remove_custom_protein_residues!()
    remove_custom_elements!()
end
