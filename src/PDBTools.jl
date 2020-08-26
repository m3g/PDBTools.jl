
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

  # Element properties
  include("./elements.jl")
  include("./atom_properties.jl")

  # For selections
  include("./select/aa_properties.jl")
  include("./select/which_natural_aminoacid.jl")
  include("./select/isacidic.jl")
  include("./select/isaliphatic.jl")
  include("./select/isaromatic.jl")
  include("./select/isbackbone.jl")
  include("./select/isbasic.jl")
  include("./select/ischarged.jl")
  include("./select/ishydrophobic.jl")
  include("./select/isneutral.jl")
  include("./select/isnonpolar.jl")
  include("./select/ispolar.jl")
  include("./select/isprotein.jl")
  include("./select/issidechain.jl")
  include("./select/iswater.jl")
  include("./select/proteins.jl")
  include("./select/select.jl")

end
