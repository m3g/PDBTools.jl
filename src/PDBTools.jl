
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
  include("./xName.jl")
  include("./xCA.jl")
  include("./xCB.jl")

  # For selections
  include("./select/isprotein.jl")
  include("./select/iswater.jl")
  include("./select/select.jl")

end
