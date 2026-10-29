#
# Methods of the functions defined in PDBToolsBase that accept selection strings. The 
# selection strings are parsed into selection functions, which are accepted by the 
# PDBToolsBase methods.
#
read_pdb(file::AbstractString, selection::AbstractString; kargs...) =
    read_pdb(file, parse_query(selection); kargs...)
read_pdb(pdbdata::IOBuffer, selection::AbstractString; kargs...) =
    read_pdb(pdbdata, parse_query(selection); kargs...)
read_mmcif(file::Union{AbstractString,IOBuffer}, selection::AbstractString; kargs...) =
    read_mmcif(file, parse_query(selection); kargs...)
write_pdb(filename::AbstractString, atoms::AbstractVector{<:Atom}, selection::AbstractString; kargs...) =
    write_pdb(filename, atoms, parse_query(selection); kargs...)
write_mmcif(filename::AbstractString, atoms::AbstractVector{<:Atom}, selection::AbstractString; kargs...) =
    write_mmcif(filename, atoms, parse_query(selection); kargs...)
read_xyz(file::Union{AbstractString,IOBuffer}, selection::AbstractString) =
    read_xyz(file, parse_query(selection))
write_xyz(filename::AbstractString, atoms::AbstractVector{<:Atom}, selection::AbstractString; kargs...) =
    write_xyz(filename, atoms, parse_query(selection); kargs...)
positions(atoms::AbstractVector{<:Atom}, selection::AbstractString) = positions(atoms, parse_query(selection))
positions(residue::Residue, selection::AbstractString) = positions(residue, parse_query(selection))
maxmin(atoms::AbstractVector{<:Atom}, selection::AbstractString) = maxmin(atoms, parse_query(selection))

@testitem "string selections" begin
    using LazyArtifacts
    using PDBTools
    ca(at) = isprotein(at) && name(at) == "CA"
    # reading
    @test read_pdb(PDBTools.test_file("structure.pdb"), "protein and name CA") == read_pdb(PDBTools.test_file("structure.pdb"), ca)
    @test length(read_pdb(PDBTools.test_file("structure.pdb"), "protein and name CA")) == 104
    @test read_pdb(PDBTools.test_file("structure.pdb"), "protein and name CA") isa Structure
    pdbdata = IOBuffer(read(PDBTools.test_file("structure.pdb"), String))
    @test length(read_pdb(pdbdata, "protein and name CA")) == 104
    @test read_mmcif(PDBTools.test_file("small.cif"), "resname HOH") == read_mmcif(PDBTools.test_file("small.cif"), at -> resname(at) == "HOH")
    @test read_mmcif(IOBuffer(read(PDBTools.test_file("small.cif"), String)), "index < 3") == read_mmcif(PDBTools.test_file("small.cif"))[1:2]
    @test read_pdb(PDBTools.test_file("small.cif"), "protein") == read_mmcif(PDBTools.test_file("small.cif"), "protein")
    @test read_xyz(PDBTools.test_file("ethanol.xyz"), "element C") == read_xyz(PDBTools.test_file("ethanol.xyz"), at -> element(at) == "C")
    @test length(read_xyz(PDBTools.test_file("ethanol.xyz"), "element C")) == 2
    @test length(read_xyz(IOBuffer(read(PDBTools.test_file("ethanol.xyz"), String)), "not element H")) == 3
    @test length(read_pdb(@view(string(PDBTools.test_file("structure.pdb"), "abc")[1:end-3]), @view("resname ALAabc"[1:end-3]))) == 72
    # writing
    pdb = read_pdb(PDBTools.test_file("small.pdb"))
    tmp = tempname()
    write_pdb(tmp, pdb, "name CA"; header=nothing)
    @test length(read_pdb(tmp)) == 3
    write_mmcif(tmp, pdb, "name CA")
    @test length(read_mmcif(tmp)) == 3
    write_xyz(tmp, pdb, "name CA")
    @test length(read_xyz(tmp)) == 3
    write_xyz(tmp, pdb, "name CA"; append=true)
    @test length(eachmodel(read_xyz(tmp))) == 2
    # coordinates
    atoms = read_pdb(PDBTools.test_file("structure.pdb"))
    @test positions(atoms, "residue = 3") == positions(atoms, at -> residue(at) == 3)
    r1 = first(eachresidue(atoms))
    @test positions(r1, "name N") == positions(r1, at -> name(at) == "N")
    @test maxmin(atoms, "protein").xmin == maxmin(atoms, isprotein).xmin
end
