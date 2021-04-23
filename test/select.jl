@testset "Selections" begin

  atoms = readPDB("./structure.pdb")
  
  @test length(select(atoms,"name CA")) == 104
  sel = select(atoms,"index = 13") 
  @test length(sel) == 1
  @test sel[1].index == 13
  @test sel[1].index_pdb == 13

  @test length(select(atoms,"index > 1 and index < 13")) == 11

  @test length(select(atoms,"protein")) == 1463

  @test length(select(atoms,"water")) == 58014
   
  @test length(select(atoms,"resname GLY")) == 84

  @test length(select(atoms,"segname PROT")) == 1463

  @test length(select(atoms,"residue = 2")) == 11

  @test length(select(atoms,"neutral")) == 1182

  @test length(select(atoms,"charged")) == 281

  @test length(select(atoms,"sidechain")) == 854

  @test length(select(atoms,"acidic")) == 162

  @test length(select(atoms,"basic")) == 68

  @test length(select(atoms,"hydrophobic")) == 327

  @test length(select(atoms,"hydrophobic")) == 327

  @test length(select(atoms,"aliphatic")) == 379

  @test length(select(atoms,"aromatic")) == 344

  @test length(select(atoms,"polar")) == 880

  @test length(select(atoms,"nonpolar")) == 583

  @test maxmin(atoms,"chain A").xlength ≈ [83.083, 83.028, 82.7] 
 
  # Test editing a field
  atoms[1].index = 0
  @test atoms[1].index == 0

  # Test residue iterator
  someresidues = select(atoms,"residue < 15")
  n = 0
  m = 0.
  for res in eachresidue(someresidues)
    if name(res) == "SER"
      n += 1
      for atom in res
        m += mass(atom)
      end
    end
  end
  @test n == 4
  @test m ≈ 348.31340000000006

  # Residue properties (discontinous set)
  lessresidues = select(someresidues,"residue < 3 or residue > 12")
  @test residue.(eachresidue(lessresidues)) = [1, 2, 13, 14]
  @test resnum.(eachresidue(lessresidues)) = [1, 2, 13, 14]
  @test name.(eachresidue(lessresidues)) = ["ALA", "CYS", "SER", "SER"]
  @test resname.(eachresidue(lessresidues)) = ["ALA", "CYS", "SER", "SER"]
  @test chain.(eachresidue(lessresidues)) = ["A", "A", "A", "A"]
  @test model.(eachresidue(lessresidues)) = [1, 1, 1, 1]

end
