module PDBTools

  using Printf
  using InteractiveUtils

  export readPDB, writePDB, getseq, wget, edit!, oneletter 
  export Atom, printatom, name
  export Residue, eachresidue, resname, residue, resnum, 
         chain, model, segname
  export coor, maxmin
  export atomic_number, element, mass, element_name

  #
  # Data structures
  #
  include("./empty_struct.jl")
  include("./Atom.jl")
  include("./Indexes_mmCIF_fields.jl")
  include("./Residue.jl")

  #
  # Input and output functions
  #
  include("./all.jl")
  include("./oneletter.jl")
  include("./check_mmCIF.jl")
  include("./parse_int.jl")
  include("./read_atom.jl")
  include("./write_atom.jl")
  include("./readPDB.jl")
  include("./edit.jl")
  include("./writePDB.jl")
  include("./getseq.jl")
  include("./coor.jl")
  include("./maxmin.jl")
  include("./same_residue.jl")
  include("./wget.jl")

  # Element properties
  include("./elements.jl")
  include("./atom_properties.jl")

  # For selections
  include("./select/select_includes.jl")
  include("./alternate_conformation.jl")

end
