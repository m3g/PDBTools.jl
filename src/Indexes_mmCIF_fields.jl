#
# Structure to read mmCIF field data indexes
#

mutable struct Indexes_mmCIF_fields
  index :: Int64
  name :: Int64
  resname:: Int64
  chain :: Int64
  resnum :: Int64
  x :: Int64
  y :: Int64
  z :: Int64
  b :: Int64
  occup :: Int64
end
Indexes_mmCIF_fields() = empty_struct(Indexes_mmCIF_fields)

