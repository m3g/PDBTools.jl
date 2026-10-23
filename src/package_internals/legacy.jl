#
# Some definitions to keep legacy compatibility
#
export readPDB, writePDB
const readPDB = read_pdb
writePDB(atoms::AbstractVector{<:Atom}, filename::AbstractString, args...; kargs...) =
    write_pdb(filename, atoms, args...; kargs...)

export atomic_sasa
const atomic_sasa = sasa_particles

export coor
coor(args...) = positions(args...)
coor(at::Atom) = position(at)