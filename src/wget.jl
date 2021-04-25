#
# Get PDB file from the database
#
function wget(pdb_id::String, selection::String )  
  query = parse_query(selection)
  return wget(pdb_id, only = atom -> apply_query(query,atom) )
end

function wget(pdb_id::String; only = all)
  file = download("https://files.rcsb.org/download/$(pdb_id).pdb")
  return readPDB(file,only=only)
end

