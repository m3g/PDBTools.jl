function edit!(atoms :: Vector{Atom})
  tmp_file_name = tempname()
  writePDB(atoms,tmp_file_name)
  InteractiveUtils.edit(tmp_file_name)
  read_again = readPDB(tmp_file_name)
  resize!(atoms,length(read_again))
  @. atoms = read_again
  rm(tmp_file_name)
end
