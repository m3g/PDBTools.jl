#
# Structure that contains the atom properties. It is mutable, so it can be edited.
#

mutable struct Atom
  index :: Int64 # The sequential index of the atoms in the file
  index_pdb :: Int64 # The index as written in the PDB file (might be anything)
  name :: String
  resname :: String
  chain :: String
  resnum :: Int64 # Number of residue as written in PDB file
  residue :: Int64 # Sequential residue (molecule) number in file
  x :: Float64
  y :: Float64
  z :: Float64
  b :: Float64
  occup :: Float64
  model :: Int64
  segname :: String # Segment name (cols 73:76)
end

Atom() = empty_struct(Atom)

export Atom

print_atom_title() = 
  @printf("%8s %4s %7s %5s %8s %8s %8s %8s %8s %5s %5s %5s %7s %9s\n",
          "index","name","resname","chain","resnum","residue","x","y","z","b","occup","model","segname","index_pdb") 
print_atom_line(atom :: Atom) =
  @printf("%8i %4s %7s %5s %8i %8i %8.3f %8.3f %8.3f %5.2f %5.2f %5i %7s %9i\n",
           atom.index, atom.name, atom.resname, atom.chain, atom.resnum, atom.residue, 
           atom.x, atom.y, atom.z, atom.b, atom.occup, atom.model, atom.segname, atom.index_pdb)

#
# Print a formatted list of atoms
#
function print_short_atom_list(atoms::AbstractVector{Atom})
  print_atom_title()
  for i in 1:min(length(atoms),3)
    atom = atoms[i]
    print_atom_line(atom)
  end
  if length(atoms) > 7
    @printf("%57s\n","⋮ ")
  end
  for i in max(4,length(atoms)-2):length(atoms)
    atom = atoms[i]
    print_atom_line(atom)
  end
end

function Base.show( io :: IO, atom :: Atom)
  println("   $(typeof(atom)) with fields:")
  print_atom_title()
  print_atom_line(atom)
end

function Base.show( io :: IO,::MIME"text/plain", atoms :: AbstractVector{Atom} )
  println("   Array{Atoms,1} with $(length(atoms)) atoms with fields:")
  print_short_atom_list(atoms)
end
