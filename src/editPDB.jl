#
# Reads PDB file atom data into a mutable array, such that the data can be edited
#


function editPDB(file :: String, selection :: String)
  query = parse_query(selection)
  return editPDB(file, only = atom -> apply_query(query,atom))
end

function editPDB(file :: String; only = all)
  pdb = readPDB(file, only = only)
  mutpdb = Vector{MutableAtom}(undef,length(pdb))
  @. mutpdb = MutableAtom(pdb)
end
