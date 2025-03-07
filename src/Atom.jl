"""
    Atom::DataType

Structure that contains the atom properties. It is mutable, so its fields can be modified.

Fields:

    mutable struct Atom{CustomType}
        index::Int32 # The sequential index of the atoms in the file
        index_pdb::Int32 # The index as written in the PDB file (might be anything)
        name::String7 # Atom name
        resname::String7 # Residue name
        chain::String3 # Chain identifier
        resnum::Int32 # Number of residue as written in PDB file
        residue::Int32 # Sequential residue (molecule) number in file
        x::Float32 # x coordinate
        y::Float32 # y coordinate
        z::Float32 # z coordinate
        beta::Float32 # temperature factor
        occup::Float32 # occupancy
        model::Int32 # model number
        segname::String7 # Segment name (cols 73:76)
        pdb_element::String3 # Element symbol string (cols 77:78)
        charge::Float32 # Charge (cols: 79:80)
        custom::CustomType # Custom fields
    end

### Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB)
   Vector{Atom{Nothing}} with 35 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
       2 1HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
⋮
      34    C     ASP     A        3        3   -2.626  -10.480   -7.749  1.00  0.00     1    PROT        34
      35    O     ASP     A        3        3   -1.940  -10.014   -8.658  1.00  0.00     1    PROT        35

julia> resname(atoms[1])
"ALA"

julia> chain(atoms[1])
"A"

julia> element(atoms[1])
"N"

julia> mass(atoms[1])
14.0067

julia> position(atoms[1])
3-element StaticArraysCore.SVector{3, Float32} with indices SOneTo(3):
  -9.229
 -14.861
  -5.481
```

The `pdb_element` and `charge` fields, which are frequently left empty in PDB files, are not printed. 
The direct access to the fields is considered part of the interface.

Custom fields can be set on `Atom` construction with the `custom` keyword argument. The Atom structure
will then be parameterized with the type of `custom`. 

### Example

```jldoctest
julia> using PDBTools

julia> atom = Atom(index = 0; custom=Dict(:c => "c", :index => 1));

julia> typeof(atom)
Atom{Dict{Symbol, Any}}

julia> atom.custom
Dict{Symbol, Any} with 2 entries:
  :index => 1
  :c     => "c"

julia> atom.custom[:c]
"c"
```

"""
mutable struct Atom{CustomType}
    index::Int32 # The sequential index of the atoms in the file
    index_pdb::Int32 # The index as written in the PDB file (might be anything)
    name::String7
    resname::String7
    chain::String3
    resnum::Int32 # Number of residue as written in PDB file
    residue::Int32 # Sequential residue (molecule) number in file
    x::Float32
    y::Float32
    z::Float32
    beta::Float32
    occup::Float32
    model::Int32
    segname::String7 # Segment name (cols 73:76)
    pdb_element::String3
    charge::Float32
    custom::CustomType
end

#
# Default constructor
#
function Atom(;custom::CustomType=nothing, kargs...) where {CustomType}
    atom = Atom{CustomType}(0,0,"X","XXX","X",0,0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0,"","X",0.0f0,custom)
    kargs_values = values(kargs)
    kargs_keys = keys(kargs_values)
    ntuple(length(kargs_values)) do i 
        @inline 
        f = kargs_keys[i]
        v = kargs_values[f]
        setfield!(atom, f, fieldtype(typeof(atom), f)(v))
    end
    return atom
end

#
# Constructor without custom::Nothing
#
Atom{Nothing}(;kargs...) = Atom(;custom=nothing, kargs...)

@testitem "Atom constructors" begin
    atref = Atom{Nothing}(0,0,"X","XXX","X",0,0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0,"","X",0.0f0,nothing)
    at = Atom()
    @test Base.summarysize(at) == 80
    @test all((getfield(at, f) == getfield(atref, f) for f in fieldnames(Atom)))
    at1 = Atom{Nothing}(;index=1, name="CA")
    at2 = Atom(;custom=nothing, index=1, name="CA")
    @test all((getfield(at1, f) == getfield(at2, f) for f in fieldnames(Atom)))
    @test (@allocations at = Atom()) <= 1 
    @test (@allocations at = Atom(; index=1, residue=1, name="CA")) <= 1
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
charge(atom::Atom) = atom.charge
pdb_element(atom::Atom) = atom.pdb_element

@testitem "Atom default fields" begin
    using PDBTools
    atoms = read_pdb(PDBTools.TESTPDB, "protein and residue 2")
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
    @test charge(atom) == 0.0f0
end

function Base.copy(atom::Atom{CustomType}) where {CustomType}
    if ismutabletype(CustomType)
        throw(ArgumentError("""\n
            The Atom object contains a mutable custom field of type $CustomType. To create an independent copy of it, use the `deepcopy` function. 
        
        """))
    else
        Atom{CustomType}((getfield(atom, f) for f in fieldnames(Atom))...)
    end
end

@testitem "copy atom" begin
    using PDBTools
    at1 = Atom()
    at2 = copy(at1)
    @test all(getfield(at1,f) == getfield(at2,f) for f in fieldnames(Atom))
    at1 = Atom(custom=1.0)
    at2 = copy(at1)
    @test all(getfield(at1,f) == getfield(at2,f) for f in fieldnames(Atom))
    at1 = Atom(custom=Int[])
    @test_throws ArgumentError copy(at1)
    at2 = deepcopy(at1)
    @test all(getfield(at1,f) == getfield(at2,f) for f in fieldnames(Atom))
end

"""
    add_custom_field(atom::Atom, value)

Adds a custom field to an `Atom` structure, returning a new `Atom` structure with the custom field added.
The returning Atom structure is parameterized with the type of `value`.

"""
function add_custom_field(atom::Atom, value)
    new_atom = Atom(;custom=value)
    for field in fieldnames(Atom)
        field == :custom && continue
        setproperty!(new_atom, field, getproperty(atom, field))
    end
    return new_atom
end

@testitem "Atom custom fields" begin
    atom = Atom(; custom=Dict(:a => 1, :b => "b", :index => 1))
    @test atom.index == 0
    @test atom.custom[:a] == 1
    @test atom.custom[:b] == "b" 
    at = Atom()
    at2 = add_custom_field(at, Dict(:a => 1, :b => "b", :index => 1))
    @test atom.index == 0
    @test atom.custom[:a] == 1
    @test atom.custom[:b] == "b" 
end

#
# Compatibility with AtomsBase interface
#
"""
    atomic_symbol(atom::Atom)

Returns the atomic symbol of an atom given the `Atom` structure.

"""
atomic_symbol(atom::Atom) = element_symbol(atom)

"""
    atomic_mass(atom::Atom)

Returns the atomic mass of an atom given the `Atom` structure.

"""
atomic_mass(atom::Atom) = mass(atom)

"""
    position(atom::Atom)

Returns the position of an atom given the `Atom` structure.

"""
Base.position(atom::Atom) = SVector(atom.x, atom.y, atom.z)

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
atom_line(atom::Atom; indent=0) = repeat(' ', indent)*@sprintf(
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

@testitem "atom_line" begin
    using PDBTools
    atoms = read_pdb(PDBTools.SMALLPDB, "protein and index 1")
    @test PDBTools.atom_line(atoms[1]) == 
        "       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1"
    buff = IOBuffer()
    printatom(buff, atoms[1])
    @test length(split(String(take!(buff)))) == 28
    printatom(buff, atoms[1]; compact=true) 
    @test String(take!(buff)) == "Atom(1N-ALA1A)"
    # just test reaching this line
    @test printatom(atoms[1]; compact=true) === nothing
end

"""
    printatom(atom::Atom)
    printatom(io::IO, atom::Atom)

Prints an `Atom` structure in a human-readable format, with a title line. By default the output is printed to `stdout`,
and the `io` argument can be used to specify a different output stream.

### Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.TESTPDB, "protein and residue 2");

julia> printatom(atoms[1])
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
      13    N     CYS     A        2        2   -6.351  -14.461   -5.695  1.00  0.00     1    PROT        13

julia> atoms[1] # default show method
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
      13    N     CYS     A        2        2   -6.351  -14.461   -5.695  1.00  0.00     1    PROT        13
```

"""
function printatom(io::IO, at::Atom; compact=false, title=!compact, newline=false, indent=0)
    title && println(io, atom_title)
    ln = newline ? '\n' : ""
    if !compact
        print(io, atom_line(at; indent), ln)
    else
        print(io, repeat(' ', indent)*"Atom($(index(at))$(name(at))-$(resname(at))$(resnum(at))$(chain(at)))", ln)
    end
end
printatom(atom::Atom; kargs...) = printatom(stdout, atom; kargs...)

function Base.show(io::IO, at::Atom)
    if get(io, :compact, false)::Bool
        printatom(io, at; compact=true)
    else
        printatom(io, at; compact=false)
    end
end

function Base.show(io::IO, ats::AbstractVector{<:Atom}; compact=nothing, indent=4, type=true, title=true)
    lines, cols = displaysize(io)
    haskey(ENV, "LINES") && (lines = parse(Int,ENV["LINES"]))
    haskey(ENV, "COLUMNS") && (cols = parse(Int,ENV["COLUMNS"]))
    natprint = min(lines-5, length(ats))
    io_compact = get(io, :compact, false)::Bool
    if !io_compact && cols >= 115 && !(compact == true) && lines > 4 
        type && println(io, "   $(typeof(ats)) with $(length(ats)) atoms with fields:")
        title && println(io, atom_title)
        idot = div(natprint,2) + 1
        dots = length(ats) > natprint
        for i in 1:natprint-1
            if dots && i == idot
                println(io, "⋮")
            else
                if i <= idot
                    printatom(io, ats[i]; compact=false, title=false, newline=true)
                else
                    printatom(io, ats[end-natprint+i]; compact=false, title=false, newline=true)
                end
            end
        end
        printatom(io, ats[end]; compact=false, title=false, newline=false)
    else # compact vector printing
        maxatcols = max(1,min(length(ats), div(cols, 25)))
        print(io,repeat(' ', indent)*"[ ")
        for i in 1:maxatcols - 1
            printatom(io, ats[i]; compact=true, title=false, newline=false)
            print(io, ", ")
        end
        if length(ats) <= maxatcols
            printatom(io, ats[maxatcols]; compact=true, title=false, newline=false)
            print(io, " ]")
        else
            printatom(io, ats[maxatcols]; compact=true, title=false, newline=false)
            print(io, "…")
        end
    end
end
Base.show(io::IO, ::MIME"text/plain", ats::AbstractVector{<:Atom}) = show(io, ats)

function Base.show(io::IO, ::MIME"text/plain", vecat::AbstractVector{<:AbstractVector{<:Atom}})
    println(io, typeof(vecat), "[ ")
    for v in vecat
        show(io, v; compact=true, indent=4)
        println(io)
    end
    print(io, "]")
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
    pdb = read_pdb(PDBTools.TESTPDB)
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
    return (atom1.resnum == atom2.resnum) & 
           (atom1.model == atom2.model) &
           (atom1.chain == atom2.chain) &
           (atom1.resname == atom2.resname) &
           (atom1.segname == atom2.segname)
end

@testitem "same_residue" begin
    pdb = read_pdb(PDBTools.TESTPDB, "protein")
    import PDBTools: same_residue
    @test same_residue(pdb[1], pdb[2])
    @test !same_residue(pdb[1], pdb[50])
    at1 = Atom()
    at2 = Atom()
    @test same_residue(at1, at2)
    at2.segname = "B"
    @test !same_residue(at1, at2)
end

#
# Atom elemental properties
#
"""
    element(atom::Atom)

Returns the element symbol, as a string, of an atom given the `Atom` structure.
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
    # if there is match, just return the name
    element_name = name(atom)
    iel = searchsortedfirst(element_names, element_name)
    if iel <= length(element_names) && element_name == element_names[iel]
        return elements[element_name].symbol_string 
    end
    # Now try to inferr from the atom name
    element_name = name(atom)
    if isempty(element_name) || element_name == "X"
        return nothing
    end
    # Check if the first character is number
    i0 = 1 + isdigit(first(element_name))
    imatch = searchsortedfirst(element_names, @view(element_name[i0:i0]); by=first)
    lmatch = searchsortedlast(element_names, @view(element_name[i0:i0]); by=first)
    for iel in imatch:lmatch
        el = element_names[iel]
        if lastindex(element_name) >= i0+length(el)-1 && el == @view(element_name[i0:i0+length(el)-1])
            return el == "X" ? nothing : el
        end
    end
    return nothing
end

@testitem "get element" begin
    using PDBTools
    atoms = read_pdb(PDBTools.TESTPDB, "protein and residue 2")
    @test element(atoms[1]) == "N"
    @test element(Atom()) === "X"
    @test element(Atom(pdb_element="")) === "X"
    @test element(Atom(pdb_element="N")) == "N"
    @test element(Atom(name = "N", pdb_element="X")) == "N"
    @test element(Atom(name = "X", pdb_element="A")) == "A" 
    @test element(Atom(name = "N", pdb_element="A")) == "A" 
    @test element(Atom(name = "A")) === nothing
    @test element(Atom(name = " ")) === nothing
    atom = Atom(name="", pdb_element="")
    @test element(atom) === nothing
end

#
# Auxiliary function to retrive another property for matching elements
#
get_element_property(at::Atom, property::Symbol) = get_element_property(at, Val(property))
function get_element_property(at::Atom, ::Val{property}) where {property}
    el = element(at)
    if isnothing(el) || !haskey(elements, el)
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
    element_symbol_string(atom::Atom)

Returns a string with the symbol of the element, given the `Atom` structure.

### Example

```jldoctest
julia> using PDBTools 

julia> at = Atom(name="NT3");

julia> element_symbol_string(at)
"N"
```

"""
element_symbol_string(at::Atom) = get_element_property(at, :symbol_string)

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
   mass = if element(at) == "X"
        nothing
   else
       get_element_property(at, Val(:mass))
   end
   return mass
end
function mass(atoms::AbstractVector{<:Atom}) 
    totmass = 0.0
    for at in atoms
        if isnothing(mass(at))
            throw(ArgumentError("Atom $(name(at)) does not have a mass defined"))
        end
        totmass += mass(at)
    end
    return totmass
end

#function Base.show(io::IO, atoms::AbstractVector{Atom})
#    println(io, " Structure file with ", length(atoms), " atoms. ")
#end

@testitem "fetch atomic element properties" setup=[AllocTest] begin
    using PDBTools
    using BenchmarkTools
    using .AllocTest: Allocs

    at = Atom(name="NT3")
    @test atomic_number(at) == 7
    @test element(at) == "N"
    @test element_name(at) == "Nitrogen"
    @test mass(at) == 14.0067
    @test mass([at, at]) == 28.0134
    atoms = read_pdb(PDBTools.TESTPDB, "protein")
    @test mass(atoms) ≈ 11079.704440000156
    @test mass(Atom(name="X")) === nothing
    @test mass(Atom(name=" ")) === nothing
    @test mass(Atom(name="A")) === nothing
    @test_throws ArgumentError mass([Atom(name="X")])
    @test element(Atom(name="CAL")) == "Ca"
    @test atomic_number(Atom(name="CAL")) == 20
    @test element(Atom(name="CAL", pdb_element="CA")) == "CA"
    @test atomic_number(Atom(name="CAL", pdb_element="CA")) === nothing
    a = @benchmark sum($mass, $atoms) samples=1 evals=1
    @test a.allocs == Allocs(0)
end

@testitem "AtomsBase interface" begin
    using PDBTools
    import StaticArrays
    at = Atom(name="NT3")
    @test atomic_number(at) == 7
    @test atomic_symbol(at) == :N
    @test atomic_mass(at) ≈ 14.0067
    @test position(at) ≈ StaticArrays.SVector(0.0, 0.0, 0.0)
end

@testitem "atom - show" begin
    using PDBTools
    using ShowMethodTesting
    ENV["LINES"] = 10
    ENV["COLUMNS"] = 120
    at = Atom(;segname="X")
    @test parse_show(at) ≈ """
       index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       0    X     XXX     X        0        0    0.000    0.000    0.000  0.00  0.00     0       X         0
    """
    # regex keeps only the first three lines
    @test parse_show([at for _ in 1:50]; repl=Dict(r"^((?:[^\n]*\n){3}).*"s => s"\1", r"PDBTools." => "")) ≈ """
       Vector{Atom{Nothing}} with 50 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       0    X     XXX     X        0        0    0.000    0.000    0.000  0.00  0.00     0       X         0
    """
    @test parse_show(Dict(1 => [at, at, at]); repl=Dict("PDBTools." => "")) ≈ """
    Dict{Int64, Vector{Atom{Nothing}}} with 1 entry:
       1 =>     [ Atom(0X-XXX0X), Atom(0X-XXX0X), Atom(0X-XXX0X) ]
    """
    @test parse_show([[at, at], [at, at, at]]; repl=Dict("PDBTools." => "")) ≈ """
    Vector{Vector{Atom{Nothing}}}[ 
        [ Atom(0X-XXX0X), Atom(0X-XXX0X) ]
        [ Atom(0X-XXX0X), Atom(0X-XXX0X), Atom(0X-XXX0X) ]
    ]
    """
end

