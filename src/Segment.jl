"""
    Segment(atoms::AbstractVector{<:Atom}, range::UnitRange{Int})

Segment data structure. It contains two fields: `atoms` which is a vector of
`Atom` elements, and `range`, which indicates which atoms of the `atoms` vector
compose the segment.

The Segment structure carries the properties of the segment 
it contains, but it does not copy the original vector of atoms, only the segment
meta data and the reference to the original vector.

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
struct Segment{T<:Atom,Vec<:AbstractVector{T}} <: AbstractVector{T}
    atoms::Vec
    range::UnitRange{Int}
    name::String7
end
name(segment::Segment) = segment.name
segname(segment::Segment) = segment.name
mass(segment::Segment) = mass(@view segment.atoms[segment.range])

function Segment(atoms::AbstractVector{<:Atom}, range::UnitRange{Int})
    i = range[begin]
    # Check if the range effectively corresponds to a single residue (unsafe check)
    for j = range[begin]+1:range[end]
        if atoms[j].segname != atoms[i].segname
            throw(ArgumentError("""\n 
                Range $range does not correspond to a single segment.

            """))
        end
    end
    Segment(
        atoms,
        range,
        atoms[i].segname,
    )
end
Segment(atoms::AbstractVector{<:Atom}) = Segment(atoms, 1:length(atoms))

function Base.getindex(segment::Segment, i::Int)
    i >= 0 || throw(ArgumentError("Index must be in 1:$(length(segment))"))
    (i <= length(segment)) || throw(ArgumentError("Segment has $(length(segment)) atoms, tried to fetch index $i."))
    i = first(segment.range) + i - 1
    segment.atoms[i]
end

#
# Structure and function to define the eachsegment iterator
#
struct EachSegment{T<:AbstractVector{<:Atom}}
    atoms::T
end

"""
    eachsegment(atoms::AbstractVector{<:Atom})

Iterator for the segments of a selection.

### Example

```julia-repl
julia> atoms = wget("1LBD");

julia> length(eachresidue(atoms))
238

julia> for res in eachresidue(atoms)
         println(res)
       end
 Residue of name SER with 6 atoms.
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     SER     A      225        1   45.228   84.358   70.638 67.05  1.00     1       -         1
       2   CA     SER     A      225        1   46.080   83.165   70.327 68.73  1.00     1       -         2
       3    C     SER     A      225        1   45.257   81.872   70.236 67.90  1.00     1       -         3
       4    O     SER     A      225        1   45.823   80.796   69.974 64.85  1.00     1       -         4
       5   CB     SER     A      225        1   47.147   82.980   71.413 70.79  1.00     1       -         5
       6   OG     SER     A      225        1   46.541   82.639   72.662 73.55  1.00     1       -         6

 Residue of name ALA with 5 atoms.
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       7    N     ALA     A      226        2   43.940   81.982   70.474 67.09  1.00     1       -         7
       8   CA     ALA     A      226        2   43.020   80.825   70.455 63.69  1.00     1       -         8
       9    C     ALA     A      226        2   41.996   80.878   69.340 59.69  1.00     1       -         9
                                                      ...

```

"""
eachsegment(atoms::AbstractVector{<:Atom}) = EachSegment(atoms)

# Collect residues default constructor
Base.collect(segments::EachSegment) = collect(Segment, segments)
Base.length(segments::EachSegment) = sum(1 for segment in segments)
Base.firstindex(segments::EachSegment) = 1
Base.lastindex(segments::EachSegment) = length(segments)
function Base.getindex(::EachSegment, ::Int)
    throw(ArgumentError("""\n
        The eachsegment iterator does not support indexing. 
        Use collect(eachsegment(atoms)) to get an indexable list of segments.
    """))
end

# Array interface
Base.size(segment::Segment) = (length(segment.range),)
Base.length(segment::Segment) = length(segment.range)
Base.eltype(::Segment) = Atom


#
# Iterate, lazily, over the segments of a structure
#
same_segment(at1::Atom, at2::Atom) = at1.segname == at2.segname && at1.model == at2.model
function Base.iterate(segments::EachSegment, current_atom=firstindex(segments.atoms))
    current_atom > length(segments.atoms) && return nothing
    next_atom = current_atom + 1
    while next_atom <= length(segments.atoms) && 
          same_segment(segments.atoms[current_atom], segments.atoms[next_atom]) 
        next_atom += 1
    end
    return (Segment(segments.atoms, current_atom:next_atom-1), next_atom)
end

#
# Iterate over atoms of one residue
#
function Base.iterate(segment::Segment, current_atom=firstindex(segment))
    current_atom > lastindex(segment) && return nothing
    return (segment[current_atom], current_atom + 1)
end


@testitem "Segment iterator" begin
    using PDBTools
    atoms = read_pdb(PDBTools.DIMERPDB)
    segments = eachsegment(atoms)
    @test length(segments) == 2
    @test firstindex(segments) == 1
    @test lastindex(segments) == 2
    @test_throws ArgumentError Segment(atoms, 1904:1910) 
    @test_throws ArgumentError segments[1]
    s = collect(segments)
    @test name(s[1]) == "A"
    @test name(s[2]) == "B"
    @test s[1].range == 1:1905
    @test s[2].range == 1906:1997
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

function Base.show(io::IO, segments::EachSegment)
    print(io, " Iterator with $(length(segments)) segments.")
end

@testitem "Segment show" begin
    using PDBTools
    using ShowMethodTesting
    ENV["LINES"] = 10
    ENV["COLUMNS"] = 120
    ats = read_pdb(PDBTools.DIMERPDB)
    s = eachsegment(ats)
    @test parse_show(s) ≈ "Iterator with 2 segments."
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
