"""
    Model(atoms::AbstractVector{<:Atom}, range::UnitRange{Int})

Model data structure. It contains two fields: `atoms` which is a vector of
`Atom` elements, and `range`, which indicates which atoms of the `atoms` vector
compose the model. Models must be consecutive in the `atoms` vector, and
are identified by having the same `model` field.
 
The Model structure carries the properties of the model
it contains, but it does not copy the original vector of atoms, only the model
meta data and the reference to the original vector.

### Example

In the example below, 1RDE is PDB entry with 11 models.

```jldoctest
julia> using PDBTools

julia> ats = wget("1RDE");

julia> models = collect(eachmodel(ats))
2-element Vector{Model}[
    A-(1905 atoms))
    B-(92 atoms))
]

julia> segname.(models[1:2])
2-element Vector{InlineStrings.String7}:
 "A"
 "B"

julia> length(models[2])
92

```

"""
struct Model{T<:Atom,Vec<:AbstractVector{T}} <: AbstractStructuralElement{T}
    atoms::Vec
    range::UnitRange{Int}
    number::Int
end

# Necessary for the interface: define the _same function
_same(::Type{Model}, at1::Atom, at2::Atom) = at1.model == at2.model

# Constructors
function Model(atoms::AbstractVector{<:Atom}, range::AbstractRange{<:Integer})
    i = range[begin]
    # Check if the range effectively corresponds to a single residue (unsafe check)
    for j = range[begin]+1:range[end]
        if !(_same(Model, atoms[i], atoms[j]))
            throw(ArgumentError("""\n 
                Range $range does not correspond to a single model.

            """))
        end
    end
    Model(
        atoms,
        UnitRange{Int}(range),
        Int(atoms[i].model),
    )
end
Model(atoms::AbstractVector{<:Atom}) = Model(atoms, eachindex(atoms))

#
# Define the iterator
#
"""
    eachmodel(atoms::AbstractVector{<:Atom})

Iterator for the models of a selection.

### Example

```jldoctest
julia> using PDBTools

julia> ats = wget("1RDE");

julia> models = eachmodel(ats)
 Iterator with 11 models.

julia> for model in models
           @show length(model)
       end
length(seg) = 1905
length(seg) = 92

julia> collect(models)
2-element Vector{Segment}[ 
    A-(1905 atoms))
    B-(92 atoms))
]
```

"""
eachmodel(atoms::AbstractVector{<:Atom}) = EachStructuralElement{Model}(atoms)

# Specific getters for this type
model(model::Model) = model.number

#
# Testing interface
#
@testitem "Model iterator" begin
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
function Base.show(io::IO, mod::Model)
    compact = get(io, :compact, false)::Bool
    if compact
        print(io, "$(model(mod))-($(length(mod)) atoms))")
    else
        println(io, " Model $(model(mod)) with $(length(mod)) atoms.")
        show(IOContext(io, :type => false), @view mod.atoms[mod.range])
    end
end

function Base.show(io::IO, model_iterator::EachStructuralElement{Model})
    print(io, " Iterator with $(length(model_iterator)) models.")
end

@testitem "Model show" begin
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
