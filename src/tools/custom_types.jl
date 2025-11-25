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
    for s in ("s", "n", "c")
        add_protein_residue!("$(s)A",PDBTools.protein_residues["ALA"])
        add_protein_residue!("$(s)R",PDBTools.protein_residues["ARG"])
        add_protein_residue!("$(s)N",PDBTools.protein_residues["ASN"])
        add_protein_residue!("$(s)D",PDBTools.protein_residues["ASP"])
        add_protein_residue!("$(s)C",PDBTools.protein_residues["CYS"])
        add_protein_residue!("$(s)Q",PDBTools.protein_residues["GLN"])
        add_protein_residue!("$(s)E",PDBTools.protein_residues["GLU"])
        add_protein_residue!("$(s)G",PDBTools.protein_residues["GLY"])
        add_protein_residue!("$(s)He",PDBTools.protein_residues["HIS"])
        add_protein_residue!("$(s)I",PDBTools.protein_residues["ILE"])
        add_protein_residue!("$(s)L",PDBTools.protein_residues["LEU"])
        add_protein_residue!("$(s)K",PDBTools.protein_residues["LYS"])
        add_protein_residue!("$(s)M",PDBTools.protein_residues["MET"])
        add_protein_residue!("$(s)F",PDBTools.protein_residues["PHE"])
        add_protein_residue!("$(s)P",PDBTools.protein_residues["PRO"])
        add_protein_residue!("$(s)S",PDBTools.protein_residues["SER"])
        add_protein_residue!("$(s)T",PDBTools.protein_residues["THR"])
        add_protein_residue!("$(s)W",PDBTools.protein_residues["TRP"])
        add_protein_residue!("$(s)Y",PDBTools.protein_residues["TYR"])
        add_protein_residue!("$(s)V",PDBTools.protein_residues["VAL"])
        add_protein_residue!("$(s)X",
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
    end
#! format: on
    return nothing
end

sirah_radii = Dict{String7, Float32}([
    ["GNz", "T1"] .=> 0.55000;
    ["GOz", "T2"] .=> 0.55000;
    ["GNn", "T3"] .=> 0.40000;
    ["GOn", "T4"] .=> 0.40000;
    "GN" => 0.42000;
    "GO" => 0.42000;
    "GC" => 0.42000;
    ["Y1C",  "Y1"] .=> 0.42000;
    ["Y2Ca",  "YA"] .=>   0.42000;
    ["Y3Sm",  "Y3"] .=>   0.42000;
    ["Y4Cv",  "Y4"] .=>   0.42000;
    ["Y5Sx",  "Y5"] .=>   0.42000;
    ["Y6Cp",  "Y6"] .=>   0.42000;
    ["A1C ",  "W1"] .=>   0.35000;
    ["A1Cw",  "WC"] .=>   0.35000;
    ["A2C ",  "W2"] .=>   0.35000;
    ["A2C ",  "W2"] .=>   0.35000;
    ["A3P ",  "W3"] .=>   0.35000;
    ["A4O ",  "W4"] .=>   0.35000;
    ["A5D ",  "WD"] .=>   0.35000;
    ["A5E ",  "WE"] .=>   0.35000;
    ["A7N ",  "WN"] .=>   0.35000;
    ["A8P ",  "WH"] .=>   0.35000;
    ["P1O ",  "QO"] .=>   0.40000;
    ["P1S ",  "QS"] .=>   0.41000;
    ["P2P ",  "QP"] .=>   0.40000;
    ["P3Cn",  "QC"] .=>   0.40000;
    ["P3Cq",  "QK"] .=>   0.40000;
    ["P4O ",  "QX"] .=>   0.40000;
    ["P5N ",  "QN"] .=>   0.40000;
    ["C1Ck",  "X1"] .=>   0.40000;
    ["C2Cr",  "X2"] .=>   0.42000;
    ["C3Cr",  "X3"] .=>   0.40000;
    ["C4Cd",  "XD"] .=>   0.40000;
    ["C4Ce",  "XE"] .=>   0.40000;
    ["C5N ",  "X5"] .=>   0.45000;
    ["C6O ",  "X6"] .=>   0.45000;
    ["C7Nk",  "X7"] .=>   0.55000;
    ["C7Nm",  "XM"] .=>   0.55000;
    ["Y5Sz",  "YZ"] .=>   0.42000;
    "PL" =>   0.46327;
    "PS" =>   0.46327;
    ["ETA",  "ET"] .=>   0.35000;
    "WX" =>   0.35000;
    ["YGL",  "GL"] .=>   0.42000;
    ["P1E",  "E1"] .=>   0.41000;
    ["P2E",  "E2"] .=>   0.41000;
    ["Y2C",  "Y2"] .=>   0.42000;
    ["Y3C",  "YT"] .=>   0.42000;
    ["Y4C",  "YF"] .=>   0.42000;
    "l2" =>   0.42000;
    "k3" =>   0.42000;
    "K3" =>   0.42000;
    "K4" =>   0.42000;
    "SX" =>   0.35000;
    "SP" =>   0.35000;
    "a2" =>   0.35000;
    "2a" =>   0.35000;
    "A2" =>   0.35000;
    "A4" =>   0.35000;
    "r2" =>   0.35000;
    "R2" =>   0.35000;
    "R3" =>   0.35000;
    "R4" =>   0.35000;
    "R6" =>   0.35000;
    "R7" =>   0.35000;
    "R8" =>   0.35000;
    "PX" =>   0.46327;
    "KX" =>   0.42906;
    "KN" =>   0.33997;
    "D2" =>   0.26698;
    "D1" =>   0.32500;
    "D6" =>   0.32500;
    "J1" =>   0.32500;
    "J2" =>   0.32500;
    "M3" =>   0.32500;
    "S4" =>   0.32500;
    "S3" =>   0.32500;
    "S2" =>   0.29599;
    "M2" =>   0.29599;
    "M4" =>   0.29599;
    "J6" =>   0.29599;
    "WT" =>   0.42000;
    "WL" =>   0.65600;
    "NaW" =>   0.58000;
    "KW" =>   0.64500;
    "ClW" =>   0.68000;
    "MgX" =>   0.40000;
    "CaX" =>   0.40000;
    ["ZnX",  "ZX"] .=>   0.40000;
])

function custom_elements!(::Type{SIRAH}; elements=PDBTools.elements)
    @warn """\n
        The element masses are not the coarse-grained ones. This must be fixed in the future.

    """ _file = nothing _line = nothing
    add_element!("GN", elements["N"])#; mass=, vdw_radius=)
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
    append!(backbone_atoms, ["GN", "GC", "GO"])
    append!(not_side_chain_atoms, ["GN", "GC", "GO"])
    nothing
end

@testitem "SIRAH" begin
    using PDBTools
    pdb_file = PDBTools.SIRAHPDB
    remove_custom_protein_residues!()
    remove_custom_elements!()
    custom_protein_residues!(SIRAH)
    custom_elements!(SIRAH)
    pdb = read_pdb(pdb_file)
    @test length(isprotein.(pdb)) == 22
    @test element(pdb[1]) == "N"
    @test element(pdb[2]) == "C"
    @test element(pdb[3]) == "C"
    bb = select(pdb, "backbone")
    @test length(bb) == 15
    sc = select(pdb, "sidechain")
    @test length(sc) == 7
    remove_custom_protein_residues!()
    remove_custom_elements!()
end
