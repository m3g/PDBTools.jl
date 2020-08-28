
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
  include("./same_residue.jl")
  export readPDB, editPDB, writePDB, getseq, coor

  # Element properties
  include("./elements.jl")
  include("./atom_properties.jl")
  export atomic_number, element, mass, name

  # For selections
  include("./select/Select.jl")

end

