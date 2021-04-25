#
# Structure to read mmCIF field data indexes
#

mutable struct Indexes_mmCIF_fields
  index::Int
  name::Int
  resname::Int
  chain::Int
  resnum::Int
  x::Int
  y::Int
  z::Int
  b::Int
  occup::Int
end
Indexes_mmCIF_fields() = empty_struct(Indexes_mmCIF_fields)

