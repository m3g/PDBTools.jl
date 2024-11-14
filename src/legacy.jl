#
# Some definitions to keep legacy compatibility
#
export readPDB, writePDB
const readPDB = read_pdb
writePDB(atoms::AbstractVector{<:Atom}, filename::String, args...; kargs...) = 
    write_pdb(filename, atoms; args...; kargs...)