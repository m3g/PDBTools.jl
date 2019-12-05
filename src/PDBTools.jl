
module PDBTools

  using Printf

  #
  # Data structures
  #

  include("./empty_struct.jl")

  include("./Atom.jl")
  export Atom

  include("./Indexes_mmCIF_fields.jl")

  #
  # Input and output functions
  #

  include("./check_mmCIF.jl")
  include("./read_atom.jl")
  include("./write_atom.jl")
  include("./readPDB.jl")
  export readPDB
  export editPDB

end
