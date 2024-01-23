"""
    Atom::DataType

Structure that contains the atom properties. It is mutable, so it can be edited. 
Fields:

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
        pdb_element::String # Element symbol string (cols 77:78)
        charge::String # Charge (cols: 79:80)
        custom::Dict{Symbol, Any} # Custom fields
    end

### Example

```jldoctest
julia> using PDBTools

julia> pdb = readPDB(PDBTools.TESTPDB); # testing PDB file

julia> printatom(pdb[1])
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1

julia> resname(pdb[1])
"ALA"

julia> chain(pdb[1])
"A"

julia> element(pdb[1])
"N"

julia> mass(pdb[1])
14.0067

julia> position(pdb[1])
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  -9.229
 -14.861
  -5.481
```

The `pdb_element` and `charge` fields, which are frequently left empty in PDB files, are not printed. 
The direct access to the fields is considered part of the interface.

Custom fields can be set on `Atom` construction with the `custom` keyword argument, which receives a 
`Dict{Symbol,Any}` as parameter. They can be retrieved with the `custom_field` function or, if the custom 
field names does not overlap with an existing field, with the dot syntax. Requires PDBTools > 0.14.3.

### Example

```jldoctest
julia> using PDBTools

julia> atom = Atom(index = 0; custom=Dict(:c => "c", :index => 1));

julia> atom.c
"c"

julia> index(atom)
0

julia> custom_field(atom, :index)
1
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
    pdb_element::String = "X"
    charge::Union{Nothing,String} = nothing
    custom::Dict{Symbol,Any} = Dict{Symbol,Any}()
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
pdb_element(atom::Atom) = atom.pdb_element
charge(atom::Atom) = atom.charge
custom_field(atom::Atom, field::Symbol) = atom.custom[field]

import Base: getproperty
function getproperty(atom::Atom, field::Symbol)
    if field in fieldnames(Atom)
        getfield(atom, field)
    else
        atom.custom[field]
    end
end

@testitem "Atom default fields" begin
    using PDBTools
    atoms = readPDB(PDBTools.TESTPDB, "protein and residue 2")
    atom = atoms[1]
    @test index(atom) == 13
    @test name(atom) == "N"
    @test resname(atom) == "CYS"
    @test chain(atom) == "A"
    @test resnum(atom) == 2
    @test residue(atom) == 2
    @test all((atom.x, atom.y, atom.z) .≈ (-6.351, -14.461, -5.695))
    @test occup(atom) == 1.0
    @test beta(atom) == 0.0
    @test model(atom) == 1
    @test segname(atom) == "PROT"
    @test index_pdb(atom) == 13
end

@testitem "Atom custom fields" begin
    atom = Atom()
    @test length(atom.custom) == 0
    atom = Atom(; custom=Dict(:a => 1, :b => "b", :index => 1))
    @test atom.index == 0
    @test atom.a == 1
    @test atom.b == "b" 
    @test custom_field(atom, :a) == 1
    @test custom_field(atom, :b) == "b"
    @test custom_field(atom, :index) == 1
end

#
# Compatibility with AtomsBase interface
#
atomic_symbol(atom::Atom) = element_symbol(atom)
atomic_mass(atom::Atom) = mass(atom)
position(atom::Atom) = SVector(atom.x, atom.y, atom.z)

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

"""
    printatom(atom::Atom)

Prints an `Atom` structure in a human-readable format.

### Example

```jldoctest
julia> using PDBTools

julia> atoms = readPDB(PDBTools.TESTPDB, "protein and residue 2")
   Array{Atoms,1} with 11 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
      13    N     CYS     A        2        2   -6.351  -14.461   -5.695  1.00  0.00     1    PROT        13
      14   HN     CYS     A        2        2   -6.473  -15.272   -5.125  0.00  0.00     1    PROT        14
      15   CA     CYS     A        2        2   -5.113  -13.737   -5.466  1.00  0.00     1    PROT        15
                                                       ⋮ 
      21  HG1     CYS     A        2        2   -3.403  -16.785   -4.019  0.00  0.00     1    PROT        21
      22    C     CYS     A        2        2   -4.610  -13.207   -6.811  1.00  0.00     1    PROT        22
      23    O     CYS     A        2        2   -4.443  -13.972   -7.759  1.00  0.00     1    PROT        23

julia> printatom(atoms[1])
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
      13    N     CYS     A        2        2   -6.351  -14.461   -5.695  1.00  0.00     1    PROT        13
```

"""
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
    element(atom::Atom)

Returns the element symbol, as a string, of an atom given its name, or `Atom` structure.
If the `pdb_element` is empty or "X", the element is inferred from the atom name. 
Othwerwise, the `pdb_element` is returned.

### Example

```jldoctest
julia> using PDBTools

julia> at = Atom(name="NT3");

julia> element(at)
"N"
```

"""
function element(atom::Atom)
    # First, check if it was defined in pdb_element
    element_name = pdb_element(atom)
    if !isempty(element_name) && element_name != "X"
        return element_name
    end
    # Now try to inferr from the atom name
    element_name = name(atom)
    if isempty(element_name) || element_name == "X"
        return nothing
    end
    # if there is match, just return the name
    iel = searchsortedfirst(element_names, element_name)
    if iel <= length(element_names) && element_name == element_names[iel]
        return element_name
    end
    # Check if the first character is number
    i0 = 1 + isdigit(first(element_name))
    imatch = searchsortedfirst(element_names, element_name[i0:i0]; by=x -> x[1])
    lmatch = searchsortedlast(element_names, element_name[i0:i0]; by=x -> x[1])
    for iel in imatch:lmatch
        el = element_names[iel]
        if lastindex(element_name) >= i0+length(el)-1 && el == element_name[i0:i0+length(el)-1]
            return el == "X" ? nothing : el
        end
    end
    return nothing
end

@testitem "get element" begin
    using PDBTools
    atoms = readPDB(PDBTools.TESTPDB, "protein and residue 2")
    @test element(atoms[1]) == "N"
    @test element(Atom()) === nothing
    @test element(Atom(pdb_element="")) === nothing
    @test element(Atom(pdb_element="N")) == "N"
    @test element(Atom(name = "N", pdb_element="X")) == "N"
    @test element(Atom(name = "X", pdb_element="A")) == "A" 
    @test element(Atom(name = "N", pdb_element="A")) == "A" 
    @test element(Atom(name = "A")) === nothing
end

#
# Auxiliary function to retrive another property for matching elements
#
function get_element_property(at::Atom, property::Symbol)
    el = element(at)
    if isnothing(el)
        return nothing
    else
        return getproperty(elements[el], property)
    end
end

"""
    atomic_number(atom::Atom)

Returns the atomic number of an atom from its `Atom` structure.

### Example

```jldoctest
julia> using PDBTools

julia> at = Atom(name="NT3");

julia> atomic_number(at)
7
```

"""
atomic_number(at::Atom) = get_element_property(at, :atomic_number)

"""
    element_name(atom::Atom)

Returns the element name of an atom given its name, or `Atom` structure.

### Example

```jldoctest
julia> using PDBTools

julia> at = Atom(name="NT3");

julia> element_name(at)
"Nitrogen"
```

"""
element_name(at::Atom) = get_element_property(at, :name)


"""
    element_symbol(atom::Atom)

Returns a symbol for element name of an atom given its name, or `Atom` structure.

### Example

```jldoctest
julia> using PDBTools 

julia> at = Atom(name="NT3");

julia> element_symbol(at)
:N
```

"""
element_symbol(at::Atom) = get_element_property(at, :symbol)

"""
    mass(atom::Atom)
    mass(atoms::AbstractVector{<:Atoms})

Returns the mass of an atom given its name, or `Atom` structure, or the total mass of a vector of `Atom`s. 

If a mass is defined as a custom field in the the `Atom` structure, it is returned. Otherwise, the mass is retrieved from the
element mass as inferred from the atom name.

## Example

```jldoctest
julia> using PDBTools

julia> atoms = [ Atom(name="NT3"), Atom(name="CA") ];

julia> mass(atoms[1])
14.0067

julia> mass(atoms)
26.017699999999998
```

"""
function mass(at::Atom) 
    if haskey(at.custom, :mass)
        return at.custom[:mass]
    else
        return get_element_property(at, :mass)
    end
end
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
    at.custom[:mass] = 1.0
    @test mass(at) == 1.0
end

@testitem "AtomsBase interface" begin
    using StaticArrays
    at = Atom(name="NT3")
    @test atomic_number(at) == 7
    @test atomic_symbol(at) == :N
    @test atomic_mass(at) ≈ 14.0067
    @test position(at) ≈ SVector(0.0, 0.0, 0.0)
end

