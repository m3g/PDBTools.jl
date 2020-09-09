
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
  include("./parse_int.jl")
  include("./read_atom.jl")
  include("./write_atom.jl")
  include("./readPDB.jl")
  include("./editPDB.jl")
  include("./writePDB.jl")
  include("./getseq.jl")
  include("./coor.jl")
  include("./maxmin.jl")
  include("./same_residue.jl")
  export readPDB, editPDB, writePDB, getseq 
  export coor, maxmin

  # Element properties
  include("./elements.jl")
  include("./atom_properties.jl")
  export atomic_number, element, mass, name

  # For selections
  include("./select/select_includes.jl")
  include("./alternate_conformation.jl")

end

