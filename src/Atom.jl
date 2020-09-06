
macro SharedAtomData() 
  ex = quote
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
  esc(ex)
end

# Mutable structure used to read data only

mutable struct MutableAtom
  @SharedAtomData()
end
MutableAtom() = empty_struct(MutableAtom)

# Immutable structure used for critical code

struct Atom
  @SharedAtomData()
end

# Initialize mutable from immutable and vice-versa

Atom( atom :: MutableAtom ) = Atom([ getfield(atom,field) for field in fieldnames(Atom) ]...)
MutableAtom( atom :: Atom ) = MutableAtom([ getfield(atom,field) for field in fieldnames(Atom) ]... )

AtomType = Union{Atom,MutableAtom}
AtomVector = Union{Vector{Atom},Vector{MutableAtom}}

print_atom_title() = 
  @printf("%8s %4s %7s %5s %8s %8s %8s %8s %8s %5s %5s %5s %7s %9s\n",
          "index","name","resname","chain","resnum","residue","x","y","z","b","occup","model","segname","index_pdb") 
print_atom_line(atom :: AtomType) =
  @printf("%8i %4s %7s %5s %8i %8i %8.3f %8.3f %8.3f %5.2f %5.2f %5i %7s %9i\n",
           atom.index, atom.name, atom.resname, atom.chain, atom.resnum, atom.residue, 
           atom.x, atom.y, atom.z, atom.b, atom.occup, atom.model, atom.segname, atom.index_pdb)

function Base.show( io :: IO, atom :: AtomType )
  println("   $(typeof(atom)) with fields:")
  print_atom_title()
  print_atom_line(atom)
end

function Base.show( io :: IO,::MIME"text/plain", atoms :: AtomVector )
  println("   Array{$(typeof(atoms[1])),1} with $(length(atoms)) atoms with fields:")
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




