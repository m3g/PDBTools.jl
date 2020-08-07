
module PDBTools

  using Printf

  #
  # Data structures
  #

  include("./empty_struct.jl")
  include("./Atom.jl")
  include("./Indexes_mmCIF_fields.jl")

  #
  # Input and output functions
  #
  include("./oneletter.jl")
  include("./check_mmCIF.jl")
  include("./read_atom.jl")
  include("./write_atom.jl")
  include("./readPDB.jl")
  include("./editPDB.jl")
  include("./writePDB.jl")
  include("./getseq.jl")
  include("./coor.jl")

  # Element properties
  include("./elements.jl")
  include("./element_properties.jl")

  # For selections
  include("./select/select.jl")
  include("./select/iswater.jl")
  include("./select/proteins.jl")

end
