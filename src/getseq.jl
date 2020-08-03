#
# Functions to read the sequence of a PDB file
#

# From the vector of atoms already read

function getseq( atoms :: Union{Vector{Atom},Vector{MutableAtom}}, selection :: String)
  query = parse_query(selection)
  natoms = length(atoms)
  n = 0
  iresidue = -1
  for i in 1:natoms
    if apply_query(query,atoms[i]) && atoms[i].resnum != iresidue
      n = n + 1
      iresidue = atoms[i].resnum
    end
  end
  seq = Array{String}(undef,n,2)
  ichain = 0
  iresidue = -1
  for i in 1:natoms 
    if apply_query(query,atoms[i]) && atoms[i].resnum != iresidue
      ichain = ichain + 1
      seq[ichain,1] = atoms[i].resname
      seq[ichain,2] = oneletter(atoms[i].resname)
      iresidue = atoms[i].resnum
    end
  end
  return seq
end

# If no selection is provided, select everything that is a protein

getseq( atoms :: Union{Vector{Atom},Vector{MutableAtom}} ) = getseq(atoms,"protein")
  
# From the file name

function getseq( file :: String, selection :: String) 
  atoms = readPDB(file)
  return getseq(atoms, selection)
end
getseq(file :: String) = getseq(file,"protein")





