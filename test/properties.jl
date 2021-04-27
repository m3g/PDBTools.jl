@testset "Properties" begin

  atoms = readPDB("./structure.pdb","protein")
  
  @test mass(atoms) ≈ 11082.699070000157
  f = formula(select(atoms,"residue < 5"))
  @test f.formula == [("H", 25), ("C", 19), ("N", 4), ("O", 7), ("S", 1)]

  @test element.(select(atoms,"residue = 1")) == ["N", "H", "H", "H", "C", "H", 
                                                  "C", "H", "H", "H", "C", "O" ] 
  @test element_name.(select(atoms,"residue = 1")) == [ "Nitrogen", "Hydrogen", "Hydrogen",
                                                        "Hydrogen", "Carbon", "Hydrogen",
                                                        "Carbon", "Hydrogen", "Hydrogen",
                                                        "Hydrogen", "Carbon", "Oxygen" ]
  @test atomic_number.(select(atoms,"residue = 1")) == [ 7, 1, 1, 1, 6, 1, 6, 1,
                                                         1, 1, 6, 8 ]
  @test length(eachresidue(atoms)) == 104
  @test name(Residue(atoms,1:12)) == "ALA"
  @test Residue(atoms,1:12).range == 1:12

  m = maxmin(atoms)
  @test m.xmin ≈ [-14.18, -17.561, -15.369]
  @test m.xmax ≈ [18.694, 14.182, 15.909]
  @test m.xlength ≈ [32.873999999999995, 31.743000000000002, 31.278]
  s = select(atoms,"residue = 3")
  r = Residue(select(atoms,"residue = 3"))  
  @test coor(s) ==
          [  -4.383   -4.51    -3.903   -3.731   -4.938  -4.417  -5.543   -5.867   -5.451   -6.974   -2.626   -1.94
            -11.903  -11.263  -11.262  -12.076  -10.279  -9.552  -9.911  -10.85   -10.837  -11.289  -10.48   -10.014
             -6.849   -6.096   -8.062   -8.767   -8.612  -9.06   -7.784   -9.684  -10.863   -9.3     -7.749   -8.658 ]
  @test coor(s) == coor(r)
  @test coor(select(atoms,"residue = 3"),column_based=false)' == coor(select(atoms,"residue = 3"))
  @test PDBTools.same_residue(atoms[1],atoms[2]) == true
  @test PDBTools.same_residue(atoms[1],atoms[20]) == false

  @test residuename("glu") == "Glutamic acid"
  @test oneletter("Glu") == "E"
  @test threeletter("E") == "GLU"

end
