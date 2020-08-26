@testset "Selections" begin

  atoms = PDBTools.readPDB("./structure.pdb")
  
  @test length(PDBTools.select(atoms,"name CA")) == 104
  sel = PDBTools.select(atoms,"index = 13") 
  @test length(sel) == 1
  @test sel[1].index == 13
  @test sel[1].index_pdb == 13

  @test length(PDBTools.select(atoms,"index > 1 and index < 13")) == 11

  @test length(PDBTools.select(atoms,"protein")) == 1463

  @test length(PDBTools.select(atoms,"water")) == 58014
   
  @test length(PDBTools.select(atoms,"resname GLY")) == 84

  @test length(PDBTools.select(atoms,"segname PROT")) == 1463

  @test length(PDBTools.select(atoms,"residue = 2")) == 11

  @test length(PDBTools.select(atoms,"neutral")) == 1182

  @test length(PDBTools.select(atoms,"charged")) == 281

  @test length(PDBTools.select(atoms,"sidechain")) == 854

  @test length(PDBTools.select(atoms,"acidic")) == 162

  @test length(PDBTools.select(atoms,"basic")) == 68

  @test length(PDBTools.select(atoms,"hydrophobic")) == 327

  @test length(PDBTools.select(atoms,"hydrophobic")) == 327

  @test length(PDBTools.select(atoms,"aliphatic")) == 379

  @test length(PDBTools.select(atoms,"aromatic")) == 344

  @test length(PDBTools.select(atoms,"polar")) == 880

  @test length(PDBTools.select(atoms,"nonpolar")) == 583

end
