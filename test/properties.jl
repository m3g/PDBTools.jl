@testset "Properties" begin

    atoms = readPDB("./structure.pdb", "protein")

    @test mass(atoms) ≈ 11079.704440000156
    f = formula(select(atoms, "residue < 5"))
    @test f.formula == [("H", 25), ("C", 19), ("N", 4), ("O", 7), ("S", 1)]

    @test element.(select(atoms, "residue = 1")) ==
          ["N", "H", "H", "H", "C", "H", "C", "H", "H", "H", "C", "O"]
    @test element_name.(select(atoms, "residue = 1")) == [
        "Nitrogen",
        "Hydrogen",
        "Hydrogen",
        "Hydrogen",
        "Carbon",
        "Hydrogen",
        "Carbon",
        "Hydrogen",
        "Hydrogen",
        "Hydrogen",
        "Carbon",
        "Oxygen",
    ]
    @test atomic_number.(select(atoms, "residue = 1")) ==
          [7, 1, 1, 1, 6, 1, 6, 1, 1, 1, 6, 8]
    residues = collect(eachresidue(atoms))
    @test length(residues) == 104
    @test name(Residue(atoms, 1:12)) == "ALA"
    @test Residue(atoms, 1:12).range == 1:12

    m = maxmin(atoms)
    @test m.xmin ≈ [-14.18, -17.561, -15.369]
    @test m.xmax ≈ [18.694, 14.182, 15.909]
    @test m.xlength ≈ [32.873999999999995, 31.743000000000002, 31.278]

    s = select(atoms, "residue = 3")
    @test coor(s) == [
        -4.383 -11.903 -6.849
        -4.51 -11.263 -6.096
        -3.903 -11.262 -8.062
        -3.731 -12.076 -8.767
        -4.938 -10.279 -8.612
        -4.417 -9.552 -9.06
        -5.543 -9.911 -7.784
        -5.867 -10.85 -9.684
        -5.451 -10.837 -10.863
        -6.974 -11.289 -9.3
        -2.626 -10.48 -7.749
        -1.94 -10.014 -8.658
    ]

    s2 = select(atoms, "residue = 5")
    @test distance(s, s2) ≈ 3.6750402718881863
    x1 = coor(s)
    x2 = coor(s2)
    @test distance(x1, x2) ≈ 3.6750402718881863
    x1 = coor(s, xyz_in_cols = false)
    x2 = coor(s2, xyz_in_cols = false)
    @test distance(x1, x2, xyz_in_cols = false) ≈ 3.6750402718881863
    @test distance(residues[3], residues[5]) ≈ 3.6750402718881863

    r = Residue(select(atoms, "residue = 3"))
    @test coor(s) == coor(r)
    @test coor(select(atoms, "residue = 3")) == coor(residues[3])
    @test coor(residues[3])' == coor(residues[3], xyz_in_cols = false)
    @test PDBTools.same_residue(atoms[1], atoms[2]) == true
    @test PDBTools.same_residue(atoms[1], atoms[20]) == false

    @test residuename("glu") == "Glutamic acid"
    @test oneletter("Glu") == "E"
    @test threeletter("E") == "GLU"

    seq = "AEG"
    @test mass(Sequence(seq)) ≈ 257.2432
    seq = ["A", "E", "G"]
    @test mass(Sequence(seq)) ≈ 257.2432
    seq = ['A', 'E', 'G']
    @test mass(Sequence(seq)) ≈ 257.2432
    seq = ["ALA", "GLU", "GLY"]
    @test mass(Sequence(seq)) ≈ 257.2432
    seq = ["Alanine", "Glutamic acid", "Glycine"]
    @test mass(Sequence(seq)) ≈ 257.2432

end
