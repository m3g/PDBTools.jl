#
# Custom protein residues that should not be loaded by default
#

#
# Sirah force field protein residues and elements
#
struct SIRAH end

function sasa_particles(::Type{SIRAH}, args...; kargs...)
    sasa_particles(args...;
        atom_type=name,
        probe_radius=2.1,
        kargs...
    )
end

function custom_protein_residues!(::Type{SIRAH}; protein_residues=PDBTools.protein_residues)
    #! format: off
    for s in ("s", "n", "c")
        add_protein_residue!("$(s)A",protein_residues["ALA"])
        add_protein_residue!("$(s)R",protein_residues["ARG"])
        add_protein_residue!("$(s)N",protein_residues["ASN"])
        add_protein_residue!("$(s)D",protein_residues["ASP"])
        add_protein_residue!("$(s)C",protein_residues["CYS"])
        add_protein_residue!("$(s)Q",protein_residues["GLN"])
        add_protein_residue!("$(s)E",protein_residues["GLU"])
        add_protein_residue!("$(s)G",protein_residues["GLY"])
        add_protein_residue!("$(s)He",protein_residues["HIS"])
        add_protein_residue!("$(s)I",protein_residues["ILE"])
        add_protein_residue!("$(s)L",protein_residues["LEU"])
        add_protein_residue!("$(s)K",protein_residues["LYS"])
        add_protein_residue!("$(s)M",protein_residues["MET"])
        add_protein_residue!("$(s)F",protein_residues["PHE"])
        add_protein_residue!("$(s)P",protein_residues["PRO"])
        add_protein_residue!("$(s)S",protein_residues["SER"])
        add_protein_residue!("$(s)T",protein_residues["THR"])
        add_protein_residue!("$(s)W",protein_residues["TRP"])
        add_protein_residue!("$(s)Y",protein_residues["TYR"])
        add_protein_residue!("$(s)V",protein_residues["VAL"])
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

function custom_elements!(::Type{SIRAH}; elements=PDBTools.elements)
    add_element!("GN", elements["N"]; symbol_string="GN", mass=14.0067f0, vdw_radius=4.2) # 4.0 for N-terminal, 5.5 for charged N-terminal
    add_element!("GC", elements["C"]; symbol_string="GC", mass=12.011f0, vdw_radius=4.2)
    add_element!("BCG", elements["C"]; symbol_string="BCG", mass=12.011f0, vdw_radius=4.2)
    add_element!("GO", elements["O"]; symbol_string="GO", mass=15.9994f0, vdw_radius=4.2) # 4.0 for C-terminal; 5.5 for charged C-terminal
    add_element!("BCZ", elements["C"]; symbol_string="BCZ", mass=12.011f0, vdw_radius=4.2)
    add_element!("BNN1", elements["N"]; symbol_string="BNN1", mass=14.0067f0, vdw_radius=4.2)
    add_element!("BNN2", elements["N"]; symbol_string="BNN2", mass=14.0067f0, vdw_radius=4.2)
    add_element!("BSG", elements["S"]; symbol_string="BSG", mass=32.06f0, vdw_radius=4.1) # 4.2 for oxydized CYS
    add_element!("BCB", elements["C"]; symbol_string="BCB", mass=12.011f0, vdw_radius=4.2)
    add_element!("BOG", elements["O"]; symbol_string="BOG", mass=15.9994f0, vdw_radius=4.2) # 4.0 for THR, SER
    add_element!("BPG", elements["H"]; symbol_string="BPG", mass=1.0079f0, vdw_radius=5.5) # 4.0 for THR, SER
    add_element!("BND", elements["N"]; symbol_string="BND", mass=14.0067f0, vdw_radius=4.2)
    add_element!("BOD", elements["O"]; symbol_string="BOD", mass=15.9994f0, vdw_radius=4.2)
    add_element!("BOE1", elements["O"]; symbol_string="BOE1", mass=15.9994f0, vdw_radius=4.2)
    add_element!("BOE2", elements["O"]; symbol_string="BOE2", mass=15.9994f0, vdw_radius=4.2)
    add_element!("BCE1", elements["C"]; symbol_string="BCE1", mass=12.011f0, vdw_radius=4.2)
    add_element!("BCE2", elements["C"]; symbol_string="BCE2", mass=12.011f0, vdw_radius=4.2)
    add_element!("BCD", elements["C"]; symbol_string="BCD", mass=12.011f0, vdw_radius=4.2)
    add_element!("BSD", elements["S"]; symbol_string="BSD", mass=32.06f0, vdw_radius=3.5)
    add_element!("BNE", elements["N"]; symbol_string="BNE", mass=14.0067f0, vdw_radius=1.7)
    add_element!("BPE", elements["H"]; symbol_string="BPE", mass=1.0079f0, vdw_radius=5.5)
    add_element!("BCE", elements["C"]; symbol_string="BCE", mass=12.011f0, vdw_radius=4.2)
    add_element!("LN", elements["N"]; symbol_string="LN", mass=14.0067f0, vdw_radius=4.2)
    add_element!("LN1", elements["N"]; symbol_string="LN1", mass=14.0067f0, vdw_radius=4.2)
    add_element!("BC12", elements["C"]; symbol_string="BC12", mass=12.011f0, vdw_radius=4.2)
    add_element!("BC13", elements["C"]; symbol_string="BC13", mass=12.011f0, vdw_radius=4.2)
    add_element!("BC14", elements["C"]; symbol_string="BC14", mass=12.011f0, vdw_radius=4.2)
    add_element!("BC15", elements["C"]; symbol_string="BC15", mass=12.011f0, vdw_radius=4.2)
    add_element!("BCT1", elements["C"]; symbol_string="BCT1", mass=12.011f0, vdw_radius=4.2)
    add_element!("BC21", elements["C"]; symbol_string="BC21", mass=12.011f0, vdw_radius=4.2)
    add_element!("BC22", elements["C"]; symbol_string="BC22", mass=12.011f0, vdw_radius=4.2)
    add_element!("BC23", elements["C"]; symbol_string="BC23", mass=12.011f0, vdw_radius=4.2)
    add_element!("BC24", elements["C"]; symbol_string="BC24", mass=12.011f0, vdw_radius=4.2)
    add_element!("BCT2", elements["C"]; symbol_string="BCT2", mass=12.011f0, vdw_radius=4.2)
    add_element!("BFO1", elements["O"]; symbol_string="BFO1", mass=15.9994f0, vdw_radius=4.2)
    add_element!("BFO2", elements["O"]; symbol_string="BFO2", mass=15.9994f0, vdw_radius=4.2)
    add_element!("WN1", elements["N"]; symbol_string="WN1", mass=18.01534f0, vdw_radius=4.2)
    add_element!("WN2", elements["N"]; symbol_string="WN2", mass=18.01534f0, vdw_radius=4.2)
    add_element!("WP1", elements["H"]; symbol_string="WP1", mass=18.01534f0, vdw_radius=4.2)
    add_element!("WP2", elements["H"]; symbol_string="WP2", mass=18.01534f0, vdw_radius=4.2)
    add_element!("NaW", elements["Na"]; symbol_string="NaW", mass=140.152f0, vdw_radius=5.8)
    add_element!("KW", elements["Cl"]; symbol_string="KW", mass=147.1903f0, vdw_radius=6.45)
    add_element!("ClW", elements["Cl"]; symbol_string="ClW", mass=143.545f0, vdw_radius=6.8)
    append!(backbone_atoms, ["GN", "GC", "GO"])
    append!(not_side_chain_atoms, ["GN", "GC", "GO"])
    nothing
end

@testitem "SIRAH" begin
    using PDBTools
    remove_custom_protein_residues!()
    remove_custom_elements!()

    pdb_file = PDBTools.SIRAHPDB
    pdb = read_pdb(pdb_file)

    s0 = sasa_particles(pdb)
    @test sasa(s0) ≈ 768.41724
    custom_protein_residues!(SIRAH)
    custom_elements!(SIRAH)

    s1 = sasa_particles(pdb)
    @test sasa(s1) ≈ 768.41724
    s1 = sasa_particles(SIRAH, pdb)
    @test sasa(s1) ≈ 1535.7573
    @test sasa(s1, "backbone") ≈ 722.8112f0
    @test sasa(s1, "sidechain") ≈ 812.94617f0

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

#= 
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
=#
