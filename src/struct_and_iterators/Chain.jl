"""
    Chain

Creates a Chain data structure. Chains must be consecutive in the `atoms` vector,
and are identified by having the same `chain`, `segment`, and `model` fields.

The Chain structure carries the properties of the atoms it contains, but it does not
copy the original vector of atoms. This means that any changes made in the Chain
structure atoms, will overwrite the original vector of atoms. 

# Examples

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.CHAINSPDB);

julia> chains = collect(eachchain(ats))
4-element Vector{Chain}[
    Chain(A-48 atoms)
    Chain(B-48 atoms)
    Chain(C-48 atoms)
    Chain(A-45 atoms)
]

julia> chains[1]
 Chain A with 48 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ASP     A        1        1  133.978  119.386  -23.646  1.00  0.00     1    ASYN         1
       2   CA     ASP     A        1        1  134.755  118.916  -22.497  1.00  0.00     1    ASYN         2
⋮
      47 HD23     LEU     A        3        3  130.568  111.868  -26.242  1.00  0.00     1    ASYN        47
      48    O     LEU     A        3        3  132.066  112.711  -21.739  1.00  0.00     1    ASYN        48

julia> mass(chains[1])
353.37881000000016 

julia> model(chains[4])
2

julia> segname(chains[2])
"ASYN"
```

"""
@kwdef struct Chain{T<:Atom,Vec<:AbstractVector{T}} <: AbstractStructuralElement{T}
    atoms::Vec
    range::UnitRange{Int}
    chain::String3
    model::Int32
    segname::String7
end

# Necessary for the interface: define the _same function
_same(::Type{Chain}, at1::Atom, at2::Atom) = 
    at1.chain == at2.chain && at1.model == at2.model && at1.segname == at2.segname

# Constructors
function Chain(atoms::AbstractVector{<:Atom}, range::AbstractRange{<:Integer})
    i = first(range) 
    if any(!(_same(Chain, atoms[j], atoms[i])) for j in range)
        throw(ArgumentError("""\n 
                Range $range does not correspond to a single protein chain.
                
        """))
    end
    Chain(
        atoms = atoms,
        range = UnitRange{Int}(range),
        chain = chain(atoms[i]),
        model = model(atoms[i]),
        segname = segname(atoms[i]),
    )
end
Chain(atoms::AbstractVector{<:Atom}) = Chain(atoms, eachindex(atoms))

"""
    eachchain(atoms::AbstractVector{<:Atom})

Iterator for the chains of a selection. 

### Example

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.CHAINSPDB);

julia> eachchain(ats)
 Chain iterator with length = 4

julia> chains = collect(eachchain(ats))
4-element Vector{Chain}[
    Chain(A-48 atoms)
    Chain(B-48 atoms)
    Chain(C-48 atoms)
    Chain(A-45 atoms)
]
```

"""
eachchain(atoms::AbstractVector{<:Atom}) = EachStructuralElement{Chain}(atoms)

# Specific getters for this type
name(chain::Chain) = chain.chain
chain(chain::Chain) = chain.chain
model(chain::Chain) = chain.model
segname(chain::Chain) = chain.segname
mass(chain::Chain) = mass(@view chain.atoms[chain.range])

@testitem "Chain iterator" begin
    using PDBTools
    pdb = read_pdb(PDBTools.CHAINSPDB)
    ichains = eachchain(pdb)
    @test Chain(pdb, 1:48).range == 1:48
    @test Chain(pdb[1:48]).range == 1:48
    @test_throws ArgumentError Chain(pdb, 49:97).range 
    @test length(ichains) == 4
    @test firstindex(ichains) == 1
    @test lastindex(ichains) == 4
    @test_throws ArgumentError ichains[1]
    chains = collect(eachchain(pdb))
    @test name(chains[3]) == "C"
    @test index.(filter(at -> resname(at) == "ASP" && name(at) == "CA", chains[1])) == [2]
    @test length(findall(at -> resname(at) == "GLN", chains[1])) == 17
    @test mass(chains[1]) ≈ 353.37881000000016
    @test segname(chains[3]) == "ASYN"
    @test model(chains[4]) == 2
    @test chain(chains[4]) == "A"
    @test_throws ArgumentError chains[1][49]
end

#
# io show functions
#
function Base.show(io::IO, ch::Chain)
    compact = get(io, :compact, false)::Bool
    if compact
        print(io, "Chain($(chain(ch))-$(length(ch)) atoms)")
    else
        println(io, " Chain $(chain(ch)) with $(length(ch)) atoms.")
        show(IOContext(io, :type => false), @view ch.atoms[ch.range])
    end
end

@testitem "Chain show" begin
    using PDBTools
    using ShowMethodTesting
    ENV["LINES"] = 10
    ENV["COLUMNS"] = 120
    ats = read_pdb(PDBTools.CHAINSPDB)
    c = eachchain(ats)
    @test parse_show(c; repl=Dict("PDBTools." => "")) ≈ "Chain iterator with length = 4"
    cc = collect(c)
    @test parse_show(cc; repl=Dict("PDBTools." => "")) ≈ """
    4-element Vector{Chain}[ 
        Chain(A-48 atoms)
        Chain(B-48 atoms)
        Chain(C-48 atoms)
        Chain(A-45 atoms)
    ]
    """
    @test parse_show(cc[1]; repl=Dict(r"^((?:[^\n]*\n){3}).*"s => s"\1")) ≈ """
     Chain A with 48 atoms.
    index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ASP     A        1        1  133.978  119.386  -23.646  1.00  0.00     1    ASYN         1
       2   CA     ASP     A        1        1  134.755  118.916  -22.497  1.00  0.00     1    ASYN         2
    """
end