"""

```julia
Atom::DataType
```

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
    "beta",
    "occup",
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
    atom.beta,
    atom.occup,
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
