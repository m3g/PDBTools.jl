#
# Reads PDB file atom data
#

function readPDB(file :: String, selection :: String)
  query = parse_query(selection)
  return readPDB(file, only = atom -> apply_query(query,atom) )
end

function readPDB(file :: String; only = atom -> true)

  # Check if structure is in mmCIF format
  mmCIF, mmCIF_fields = check_mmCIF(file)

  # Read file
  pdbfile = open(file,"r")
  natoms = 0
  index = 0
  imodel = 1
  iresidue = 1
  atoms = Vector{Atom}(undef,0)
  local lastatom
  for line in eachline(pdbfile)
    if occursin("END",line)
      imodel = imodel + 1
    end
    atom = read_atom(line, mmCIF = mmCIF, mmCIF_fields = mmCIF_fields)
    if atom != nothing 
      index = index + 1
      atom.index = index
      atom.model = imodel
      if index > 1
        if ! same_residue(atom,lastatom)
          iresidue += 1
        end
      end
      atom.residue = iresidue
      if only(atom) 
        natoms = natoms + 1 
        push!(atoms,Atom(atom))
      end
      lastatom = atom
    end
  end
  close(pdbfile)
  if natoms == 0
    error(" Could not find any atom in PDB file matching the selection. ")
  end

  return atoms
end

import Base.show
function Base.show( io :: IO, atoms :: Array{Atom} )
  println(" Structure file with ", length(atoms), " atoms. ")
end


