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

end
