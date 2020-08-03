#
# Reads PDB file atom data into a mutable array, such that the data can be edited
#

function editPDB(file :: String, selection :: String)

  pdb = readPDB(file, selection)
  mutpdb = Vector{MutableAtom}(undef,length(pdb))
  @. mutpdb = MutableAtom(pdb)

end

editPDB(file :: String) = editPDB(file,"all")

