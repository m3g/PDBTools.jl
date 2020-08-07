
AtomData = quote
  index :: Int64 # The sequential index of the atoms in the file
  index_pdb :: Int64 # the index as written in the PDB file (might be anything)
  name :: String
  resname :: String
  chain :: String
  resnum :: Int64
  x :: Float64
  y :: Float64
  z :: Float64
  b :: Float64
  occup :: Float64
  model :: Int64
  atomic_number :: Int64
end

# Mutable structure used to read data only

@eval mutable struct MutableAtom
  $AtomData
end
MutableAtom() = empty_struct(MutableAtom)

# Immutable structure used for critical code

@eval struct Atom
  $AtomData
end

# Initialize mutable from immutable and vice-versa

Atom( atom :: MutableAtom ) = Atom([ getfield(atom,i) for i in 1:nfields(atom) ]...)
MutableAtom( atom :: Atom ) = MutableAtom([ getfield(atom,i) for i in 1:nfields(atom) ]... )

