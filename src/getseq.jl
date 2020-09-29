#
# Functions to read the sequence of a PDB file
#

# From the vector of atoms already read

function getseq( atoms :: Vector{Atom}, selection :: String)
  query = parse_query(selection)
  return getseq(atoms, only = atom -> apply_query(query,atom))
end

function getseq( atoms :: Vector{Atom}; only = atom -> isprotein(atom))
  natoms = length(atoms)
  n = 0
  iresidue = -1
  for i in 1:natoms
    if only(atoms[i]) && atoms[i].resnum != iresidue
      n = n + 1
      iresidue = atoms[i].resnum
    end
  end
  seq = Array{String}(undef,n,2)
  ichain = 0
  iresidue = -1
  for i in 1:natoms 
    if only(atoms[i]) && atoms[i].resnum != iresidue
      ichain = ichain + 1
      seq[ichain,1] = atoms[i].resname
      seq[ichain,2] = oneletter(atoms[i].resname)
      iresidue = atoms[i].resnum
    end
  end
  return seq
end

# From the file name

function getseq( file :: String, selection :: String) 
  atoms = readPDB(file)
  return getseq(atoms, selection)
end

getseq( file :: String; only = atom -> isprotein(atom)) = getseq(readPDB(file), only = only)

