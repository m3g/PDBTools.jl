#
# Reads PDB file atom data
#

readPDB(file :: String) = readPDB(file,"all")

function readPDB(file :: String, selection :: String)

  # Parse selection query
  query = parse_query(selection)

  # Check if structure is in mmCIF format
  mmCIF, mmCIF_fields = check_mmCIF(file)

  # Check number of atoms
  pdbfile = open(file,"r")
  natoms = 0
  imodel = 1
  for line in eachline(file)
    if occursin(line,"ENDMDL") || occursin(line,"END")
      imodel = imodel + 1
    end
    atom = read_atom(line, mmCIF = mmCIF, mmCIF_fields = mmCIF_fields, model = imodel)
    if atom != Nothing && apply_query(query,atom)
      natoms = natoms + 1
    end 
  end
  seek(pdbfile,0)
  if natoms == 0
    close(pdbfile)
    error(" Could not find any atom in PDB file. ")
  end
  atoms = Vector{Atom}(undef,natoms)

  # Read file
  iatom = 0
  imodel = 1
  for line in eachline(file)
    if occursin(line,"ENDMDL") || occursin(line,"END")
      imodel = imodel + 1
    end
    atom = read_atom(line, mmCIF = mmCIF, mmCIF_fields = mmCIF_fields, model = imodel)
    if atom != Nothing && apply_query(query,atom)
      atom.model = imodel
      iatom = iatom + 1 
      atoms[iatom] = Atom(atom)
    end
  end
  close(pdbfile)

  return atoms
end

import Base.show
function Base.show( io :: IO, atoms :: Array{Atom} )
  println(" Structure file with ", length(atoms), " atoms. ")
end





