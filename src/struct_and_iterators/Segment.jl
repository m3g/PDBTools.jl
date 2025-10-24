"""
    Segment

Segment data structure. Segments must be consecutive in the `atoms` vector, and
are identified by having the same `segname` and `model` fields.
 
The Segment structure carries the properties of the segment 
it contains, but it does not copy the original vector of atoms, only the segment
meta data and the reference to the original vector. Thus, changes
in the segment atoms will be reflected in the original vector of atoms.

### Example

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> segments = collect(eachsegment(ats))
2-element Vector{Segment}[
    A-(1905 atoms))
    B-(92 atoms))
]

julia> segname.(segments[1:2])
2-element Vector{InlineStrings.String7}:
 "A"
 "B"

julia> length(segments[2])
92

```

"""
struct Segment{T<:Atom,Vec<:AbstractVector{T}} <: AbstractStructuralElement{T}
    atoms::Vec
    range::UnitRange{Int}
    name::String7
end

# Necessary for the interface: define the same_struct_element function
same_struct_element(::Type{Segment}, at1::Atom, at2::Atom) = at1.segname == at2.segname && at1.model == at2.model

# Constructors 
function Segment(atoms::AbstractVector{<:Atom}, range::AbstractRange{<:Integer})
    _check_unique(Segment, atoms, range)
    return Segment(atoms, UnitRange{Int}(range), segname(first(atoms[range])))
end
Segment(atoms::AbstractVector{<:Atom}) = Segment(atoms, eachindex(atoms))

"""
    eachsegment(atoms::AbstractVector{<:Atom})

Iterator for the segments of a selection.

### Example

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> sit = eachsegment(ats)
 Segment iterator with length = 2

julia> for seg in sit
           @show length(seg)
       end
length(seg) = 1905
length(seg) = 92

julia> collect(sit)
2-element Vector{Segment}[ 
    A-(1905 atoms))
    B-(92 atoms))
]
```

"""
eachsegment(atoms::AbstractVector{<:Atom}) = EachStructuralElement{Segment}(atoms)

# Specific getters for this type
name(segment::Segment) = segment.name
segname(segment::Segment) = segment.name
mass(segment::Segment) = mass(@view segment.atoms[segment.range])
get_atoms(segment::Segment) = @view(segment.atoms[segment.range])

@testitem "Segment iterator" begin
    using PDBTools
    atoms = read_pdb(PDBTools.DIMERPDB)
    segments = eachsegment(atoms)
    @test length(segments) == 2
    @test firstindex(segments) == 1
    @test lastindex(segments) == 2
    @test length(last(segments)) == 92
    @test segname(last(segments)) == "B"
    @test_throws ArgumentError Segment(atoms, 1904:1910)
    @test_throws ArgumentError segments[1]
    s = collect(segments)
    @test name(s[1]) == "A"
    @test name(s[2]) == "B"
    @test s[1].range == 1:1905
    @test length(get_atoms(s[1])) == 1905
    @test s[2].range == 1906:1997
    @test name(s[1]) == "A"
    @test segname(s[1]) == "A"
    @test length(s[1]) == 1905
    @test mass(s[1]) ≈ 25222.33909999994
    @test size(s[1]) == (1905,)
    @test eltype(s[1]) == Atom
    @test sum(mass(at) for at in s[1]) ≈ 25222.33909999994
    s = Segment(atoms[1:1905])
    @test name(s) == "A"
end

#
# io show functions
#
function Base.show(io::IO, segment::Segment)
    compact = get(io, :compact, false)::Bool
    if compact
        print(io, "$(name(segment))-($(length(segment)) atoms))")
    else
        println(io, " Segment of name $(name(segment)) with $(length(segment)) atoms.")
        show(IOContext(io, :type => false), @view segment.atoms[segment.range])
    end
end

@testitem "Segment show" begin
    using PDBTools
    using ShowMethodTesting
    ENV["LINES"] = 10
    ENV["COLUMNS"] = 120
    ats = read_pdb(PDBTools.DIMERPDB)
    s = eachsegment(ats)
    @test parse_show(s; repl=Dict("PDBTools." => "")) ≈ "Segment iterator with length = 2"
    sc = collect(s)
    @test parse_show(sc; repl=Dict("PDBTools." => "")) ≈ """
    2-element Vector{Segment}[ 
    A-(1905 atoms))
    B-(92 atoms))
    """
    @test parse_show(sc[1]; repl=Dict(r"^((?:[^\n]*\n){3}).*"s => s"\1")) ≈ """
     Segment of name A with 1905 atoms.
    index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     LYS     A      211        1   52.884   24.022   35.587  1.00 53.10     1       A         1
       2   CA     LYS     A      211        1   52.916   24.598   36.993  1.00 53.10     1       A         2
    """
end
