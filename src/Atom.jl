
AtomData = quote
  index :: Int64
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
end

# Mutable structure used to read data only

@eval mutable struct ReadAtom
  $AtomData
end
ReadAtom() = empty_struct(ReadAtom)

# Immutable structure used for critical code

@eval struct Atom
  $AtomData
end

# Initialize mutable from immutable and vice-versa

Atom( atom :: ReadAtom ) = Atom([ getfield(atom,i) for i in 1:nfields(atom) ]...)
ReadAtom( atom :: Atom ) = ReadAtom([ getfield(atom,i) for i in 1:nfields(atom) ]... )

