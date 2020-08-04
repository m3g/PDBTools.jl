#
# Reads PDB file atom data
#

function readPDB(file :: String, selection :: String)

  # Parse selection query
  query = parse_query(selection)

  # Check if structure is in mmCIF format
  mmCIF, mmCIF_fields = check_mmCIF(file)

  # Read file
  pdbfile = open(file,"r")
  natoms = 0
  index = 0
  imodel = 1
  atoms = Vector{Atom}(undef,0)
  for line in eachline(pdbfile)
    if occursin("END",line)
      imodel = imodel + 1
    end
    atom = read_atom(line, mmCIF = mmCIF, mmCIF_fields = mmCIF_fields)
    if atom != nothing 
      index = index + 1
      if apply_query(query,atom)
        atom.index = index
        atom.model = imodel
        natoms = natoms + 1 
        push!(atoms,Atom(atom))
      end
    end
  end
  close(pdbfile)
  if natoms == 0
    error(" Could not find any atom in PDB file matching the selection. ")
  end

  return atoms
end

readPDB(file :: String) = readPDB(file,"all")

import Base.show
function Base.show( io :: IO, atoms :: Array{Atom} )
  println(" Structure file with ", length(atoms), " atoms. ")
end


