#
# Reads PDB file atom data
#

function readPDB(file :: String; chain :: String = "0", model :: Int64 = 0)

  # Check if structure is in mmCIF format

  mmCIF, mmCIF_fields = check_mmCIF(file)
  #println(" mmCIF format: ", mmCIF)

  # Check number of atoms

  pdbfile = open(file,"r")
  natoms = 0
  imodel = 1
  for line in eachline(file)
     if occursin(line,"ENDMDL")
       imodel = imodel + 1
     end
     atom = read_atom(line, mmCIF = mmCIF, mmCIF_fields = mmCIF_fields)
     if atom != Nothing
       if ( atom.chain == chain && atom.model == model ) ||
          ( atom.chain == chain && model == 0 ) ||
          ( chain == "0" && atom.model == model ) ||
          ( chain == "0" && model == 0 )
         natoms = natoms + 1
       end
     end 
  end
  seek(pdbfile,0)
  if natoms == 0
    error(" Could not find any atom in PDB file. ")
  end

  atoms = Vector{Atom}(undef,natoms)

  # Read file

  pdbfile = open(file,"r")
  iatom = 0
  imodel = 1
  for line in eachline(file)
     if occursin(line,"ENDMDL")
       imodel = imodel + 1
     end
     atom = read_atom(line, mmCIF = mmCIF, mmCIF_fields = mmCIF_fields)
     if atom != Nothing
       atom.model = imodel
       if ( atom.chain == chain && atom.model == model ) ||
          ( atom.chain == chain && model == 0 ) ||
          ( chain == "0" && atom.model == model ) ||
          ( chain == "0" && model == 0 )
         iatom = iatom + 1 
         atoms[iatom] = Atom(atom)
       end
     end
  end
  close(pdbfile)

  return atoms

end

import Base.show
function Base.show( io :: IO, atoms :: Array{Atom} )
  println(" Structure file with ", length(atoms), " atoms. ")
end





