"""
    Atom::DataType

Structure that contains the atom properties. It is mutable, so it can be edited. 
Fields:

```julia
mutable struct Atom
  index::Int # The sequential index of the atoms in the file
  index_pdb::Int # The index as written in the PDB file (might be anything)
  name::String # Atom name
  resname::String # Residue name
  chain::String # Chain identifier
  resnum::Int # Number of residue as written in PDB file
  residue::Int # Sequential residue (molecule) number in file
  x::Float64 # x coordinate
  y::Float64 # y coordinate
  z::Float64 # z coordinate
  beta::Float64 # temperature factor
  occup::Float64 # occupancy
  model::Int # model number
  segname::String # Segment name (cols 73:76)
end
```

### Example

```julia-repl
julia> pdb = wget("1LBD");

julia> printatom(pdb[1])
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ALA     P        1        1    2.062  -13.995   21.747  0.00  1.00     1    PROT         1

julia> pdb[1].resname
"ALA"

julia> pdb[1].chain
"P"

julia> element(pdb[1])
"N"

julia> mass(pdb[1])
14.0067

```

"""
Base.@kwdef mutable struct Atom
    index::Int = 0 # The sequential index of the atoms in the file
    index_pdb::Int = 0 # The index as written in the PDB file (might be anything)
    name::String = "X"
    resname::String = "XXX"
    chain::String = "X"
    resnum::Int = 0 # Number of residue as written in PDB file
    residue::Int = 0 # Sequential residue (molecule) number in file
    x::Float64 = 0.0
    y::Float64 = 0.0
    z::Float64 = 0.0
    beta::Float64 = 0.0
    occup::Float64 = 0.0
    model::Int = 0
    segname::String = "XXXX" # Segment name (cols 73:76)
end

index(atom::Atom) = atom.index
index_pdb(atom::Atom) = atom.index_pdb
name(atom::Atom) = atom.name
resname(atom::Atom) = atom.resname
chain(atom::Atom) = atom.chain
resnum(atom::Atom) = atom.resnum
residue(atom::Atom) = atom.residue
beta(atom::Atom) = atom.beta
occup(atom::Atom) = atom.occup
model(atom::Atom) = atom.model
segname(atom::Atom) = atom.segname

const atom_title = @sprintf(
    "%8s %4s %7s %5s %8s %8s %8s %8s %8s %5s %5s %5s %7s %9s",
    "index",
    "name",
    "resname",
    "chain",
    "resnum",
    "residue",
    "x",
    "y",
    "z",
    "occup",
    "beta",
    "model",
    "segname",
    "index_pdb"
)
atom_line(atom::Atom) = @sprintf(
    "%8i %4s %7s %5s %8i %8i %8.3f %8.3f %8.3f %5.2f %5.2f %5i %7s %9i",
    atom.index,
    atom.name,
    atom.resname,
    atom.chain,
    atom.resnum,
    atom.residue,
    atom.x,
    atom.y,
    atom.z,
    atom.occup,
    atom.beta,
    atom.model,
    atom.segname,
    atom.index_pdb
)

function printatom(atom::Atom)
    println(atom_title)
    println(atom_line(atom))
end

#
# Print a formatted list of atoms
#
function print_short_atom_list(io::IO, atoms::AbstractVector{Atom})
    println(io, atom_title)
    for i = 1:min(length(atoms), 3)
        print(io, atom_line(atoms[i]))
        i == length(atoms) || print(io, "\n")
    end
    if length(atoms) > 7
        @printf(io, "%57s\n", "⋮ ")
    end
    for i = max(4, length(atoms) - 2):length(atoms)
        print(io, atom_line(atoms[i]))
        i == length(atoms) || print(io, "\n")
    end
end

function Base.show(io::IO, atom::Atom)
    print(io, atom_line(atom))
end

function Base.show(io::IO, ::MIME"text/plain", atoms::AbstractVector{Atom})
    println(io, "   Array{Atoms,1} with $(length(atoms)) atoms with fields:")
    print_short_atom_list(io, atoms)
end

#
# atom properties on the structure
#
export isprotein, isbackbone, issidechain
isprotein(atom::Atom) = haskey(protein_residues, atom.resname)

const backbone_atoms = ["N", "CA", "C", "O"]
isbackbone(atom::Atom; backbone_atoms=backbone_atoms) = isprotein(atom) && atom.name in backbone_atoms

const not_side_chain_atoms = ["N", "CA", "C", "O", "HN", "H", "HA", "HT1", "HT2", "HT3"]
issidechain(atom::Atom; not_side_chain_atoms=not_side_chain_atoms) = isprotein(atom) && !(atom.name in not_side_chain_atoms)

@testitem "atoms in struct" begin
    pdb = readPDB(PDBTools.TESTPDB)
    glu = select(pdb, "resname GLU")
    @test isbackbone(glu[1])
    @test !issidechain(glu[1])
    phe = select(pdb, "resname PHE")
    @test isbackbone(phe[1])
    @test !issidechain(phe[1])
    @test issidechain(phe[8])
end

#
# Function that checks if two atoms belong to the same residue
# without, of course, checking the residue counter
#
function same_residue(atom1::Atom, atom2::Atom)
    atom1.resnum != atom2.resnum && return false
    atom1.model != atom2.model && return false
    atom1.chain != atom2.chain && return false
    atom1.resname != atom2.resname && return false
    atom1.segname != atom2.segname && return false
    return true
end

@testitem "same_residue" begin
    pdb = readPDB(PDBTools.TESTPDB, "protein")
    import PDBTools: same_residue
    @test same_residue(pdb[1], pdb[2])
    @test !same_residue(pdb[1], pdb[50])
end


#
# Atom elemental properties
#
"""
    atomic_number(name::String or atom::Atom)

Returns the atomic number of an atom given its name, or `Atom` structure.

### Example

```julia-repl
julia> at = Atom(name="NT3");

julia> atomic_number(at)
7

julia> atomic_number("CA")
6

```

"""
atomic_number(atom::Atom) = atomic_number(atom.name)

"""
    element_name(name::String or atom::Atom)

Returns the element name of an atom given its name, or `Atom` structure.

### Example

```julia-repl
julia> at = Atom(name="NT3");

julia> element_name(at)
"Nitrogen"

julia> element_name("NT3")
"Nitrogen"

julia> element_name("CA")
"Carbon"

```

"""
element(atom::Atom) = element(atom.name)
element_name(atom::Atom) = element_name(atom.name)

"""
    mass(name::String or atom::Atom or Vector{Atom})

Returns the mass of an atom given its name, or `Atom` structure, or the total mass of a vector of `Atom`s. 

### Example

```julia-repl
julia> atoms = [ Atom(name="NT3"), Atom(name="CA") ];

julia> mass(atoms[1])
14.0067

julia> mass("CA")
12.011

julia> mass(atoms)
26.017699999999998

```

"""
mass(atom::Atom) = mass(atom.name)
mass(atoms::AbstractVector{Atom}) = sum(mass, atoms)

@testitem "fetch atomic element properties" begin
    at = Atom(name="NT3")
    @test atomic_number(at) == 7
    @test element(at) == "N"
    @test element_name(at) == "Nitrogen"
    @test mass(at) == 14.0067
    @test mass([at, at]) == 28.0134
    atoms = readPDB(PDBTools.TESTPDB, "protein")
    @test mass(atoms) ≈ 11079.704440000156
end