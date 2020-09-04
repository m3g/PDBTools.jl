
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
MutableAtom( atom :: Atom ) = MutableAtom([ getfield(atom,field) for field in 1:fieldnames(Atom) ]... )

AtomType = Union{Atom,MutableAtom}
AtomVector = Union{Vector{Atom},Vector{MutableAtom}}

