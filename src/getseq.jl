#
# Functions to read the sequence of a PDB file
#

# From the file name

function getseq( file :: String; chain = nothing, protein = true) 
  atoms = readPDB(file)
  return getseq(atoms, chain = chain, protein = protein)
end

# From the vector of atoms already read

function getseq( atoms :: Vector{Atom}; chain = nothing, protein = true )
  natoms = length(atoms)
  nchain = 0
  iresidue = -1
  for i in 1:natoms
    if protein && !( isprotein(atoms[i]) )
      continue
    end
    if chain == nothing || atoms[i].chain == chain
      if atoms[i].resnum != iresidue
        nchain = nchain + 1
        iresidue = atoms[i].resnum
      end
    end
  end
  seq = Array{String}(undef,nchain,2)
  ichain = 0
  iresidue = -1
  for i in 1:natoms 
    if protein && !( isprotein(atoms[i]) )
      continue
    end
    if chain == nothing || atoms[i].chain == chain
      if atoms[i].resnum != iresidue
        ichain = ichain + 1
        seq[ichain,1] = atoms[i].resname
        seq[ichain,2] = oneletter(atoms[i].resname)
        iresidue = atoms[i].resnum
      end
    end
  end
  return seq
end




