"""
    Chain(atoms::AbstractVector{<:Atom}, range::UnitRange{Int})

Creates a Chain data structure. It contains two fields: `atoms` which is a vector of
`Atom` elements, and `range`, which indicates which atoms of the `atoms` vector
compose the chain. 

The Chain structure carries the properties of the atoms it contains, but it does not
copy the original vector of atoms. This means that any changes made in the Chain
structure atoms, will overwrite the original vector of atoms. 

# Examples

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.CHAINSPDB);

julia> chain.(eachchain(ats)) 
4-element Vector{InlineStrings.String3}:
 "A"
 "B"
 "C"
 "A"

julia> chains = collect(eachchain(ats))
   Array{Chain,1} with 4 chains.

julia> chains[1]
 Chain of name A with 48 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ASP     A        1        1  133.978  119.386  -23.646  1.00  0.00     1    ASYN         1
       2   CA     ASP     A        1        1  134.755  118.916  -22.497  1.00  0.00     1    ASYN         2
       3    C     ASP     A        1        1  135.099  117.439  -22.652  1.00  0.00     1    ASYN         3
                                                       ⋮ 
      46 HD22     LEU     A        3        3  130.704  113.003  -27.586  1.00  0.00     1    ASYN        46
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
@kwdef struct Chain{T<:Atom,Vec<:AbstractVector{T}} <: AbstractVector{T}
    atoms::Vec
    range::UnitRange{Int}
    chain::String3
    model::Int32
    segname::String7
end

name(chain::Chain) = chain.chain
chain(chain::Chain) = chain.chain
model(chain::Chain) = chain.model
segname(chain::Chain) = chain.segname
mass(chain::Chain) = mass(@view chain.atoms[chain.range])

function Chain(atoms::AbstractVector{<:Atom}, range::UnitRange{<:Integer})
    i = first(range) 
    if any(atoms[j].chain != atoms[i].chain for j in range)
        throw(ArgumentError("""\n 
                Range $range does not correspond to a single protein chain."""))
    end
    Chain(
        atoms = atoms,
        range = range,
        chain = chain(atoms[i]),
        model = model(atoms[i]),
        segname = segname(atoms[i]),
    )
end

Chain(atoms::AbstractVector{<:Atom}) = Chain(atoms, 1:length(atoms))

function Base.getindex(chain::Chain, i::Integer)
    i > 0 || throw(ArgumentError("Chain index must be in 1:$(length(chain))"))
    (i <= length(chain)) || throw(ArgumentError("Chain has $(length(chain)) atoms, tried to fetch index $i."))
    # Calculate the actual index in the atoms array
    atom_index = first(chain.range) + i - 1
    return chain.atoms[atom_index]
end

#
# Structure and function to define the eachchain iterator
#
struct EachChain{T<:AbstractVector{<:Atom}}
    atoms::T
end

"""
    eachchain(atoms::AbstractVector{<:Atom})

Iterator for the chains of a selection. 

### Example

```julia-repl
julia> ats = read_pdb(PDBTools.CHAINSPDB);

julia> length(eachchain(ats))
4

julia> for chain in eachchain(ats)
           println("Chain of name \$(name(chain))")
           println(resname.(eachresidue(chain)))
       end
A
InlineStrings.String7["ASP", "GLN", "LEU"]
B
InlineStrings.String7["ASP", "GLN", "LEU"]
C
InlineStrings.String7["ASP", "GLN", "LEU"]
A
InlineStrings.String7["ASP", "GLN", "VAL"


```

"""
eachchain(atoms::AbstractVector{<:Atom}) = EachChain(atoms)

# Collect chains default constructor
Base.collect(c::EachChain) = collect(Chain, c)
Base.length(chains::EachChain) = sum(1 for chain in chains)
Base.firstindex(chains::EachChain) = 1
Base.lastindex(chains::EachChain) = length(chains)
# Base.reverse(chains::EachChain) = reverse(chains)
# #Base.last(chains::EachChain) =  first(reverse(chains))
# Base.last(c::EachChain) = @inbounds c[max(firstindex(c), end + 1 - length(c)):end]
function Base.last(chains::EachChain) 
    local last_chain
    for chain in chains 
        last_chain = chain
    end
    return last_chain
end

function Base.getindex(::EachChain, ::Integer)
    throw(ArgumentError("""\n
        The eachchain iterator does not support indexing. 
        Use collect(eachchain(atoms)) to get an indexable list of chains.

    """))
end

Base.size(chain::Chain) = (length(chain.range),)
Base.length(chain::Chain) = length(chain.range)
Base.eltype(::Chain{T}) where {T} = T

# Iterate over residues of a structure
#
function Base.iterate(chains::EachChain, current_atom=firstindex(chains.atoms))
    current_atom > length(chains.atoms) && return nothing
    next_atom = current_atom + 1
    while next_atom <= length(chains.atoms) && 
          same_chain(chains.atoms[current_atom], chains.atoms[next_atom]) 
        next_atom += 1
    end
    return (Chain(chains.atoms, current_atom:next_atom-1), next_atom)
end

# Iterate over atoms of one residue
#
function Base.iterate(chain::Chain, current_atom=nothing)
    first_atom = index(first(chain))
    last_atom = index(last(chain))
    if isnothing(current_atom)
        current_atom = first_atom
    elseif current_atom > last_atom
        return nothing
    end
    return (chain[current_atom - first_atom + 1], current_atom + 1)
end

function same_chain(atom1::Atom, atom2::Atom)
    atom1.chain == atom2.chain &&
    atom1.model == atom2.model &&
    atom1.segname == atom2.segname
end

# io show functions
#
function Base.show(io::IO, ::MIME"text/plain", chain::Chain)
    println(io, " Chain of name $(chain.chain) with $(length(chain)) atoms.")
    print_short_atom_list(io, @view chain.atoms[chain.range])
end

function Base.show(io::IO, chains::EachChain)
    print(io, " Iterator with $(length(chains)) chains.")
end

function Base.show(io::IO, ::MIME"text/plain", chains::AbstractVector{Chain})
    print(io, "   $(typeof(chains)) with $(length(chains)) chains.")
end

@testitem "Chain iterator" begin
    pdb = read_pdb(PDBTools.CHAINSPDB)
    ichains = eachchain(pdb)
    @test Chain(pdb, 1:48).range == 1:48
    @test_throws ArgumentError Chain(pdb, 49:97).range 
    @test length(ichains) == 4
    @test firstindex(ichains) == 1
    @test lastindex(ichains) == 4
    @test last(ichains).model == 2
    @test_throws ArgumentError ichains[1]
    chains = collect(eachchain(pdb))
    @test name(chains[3]) == "C"
    @test index.(filter(at -> resname(at) == "ASP" && name(at) == "CA", chains[1])) == [2]
    @test length(findall(at -> resname(at) == "GLN", chains[1])) == 17
    @test mass(chains[1]) == 353.37881000000016
    @test segname(chains[3]) == "ASYN"
    @test model(chains[4]) == 2
    @test chain(chains[4]) == "A"
    @test_throws ArgumentError chains[1][49]
    buff = IOBuffer()
    show(buff, MIME"text/plain"(), chains[1])
    @test length(split(String(take!(buff)))) == 14*7 + 8
    show(buff, MIME"text/plain"(), ichains)
    @test String(take!(buff)) == " Iterator with 4 chains."
    show(buff, MIME"text/plain"(), chains)
    @test String(take!(buff)) == "   Vector{PDBTools.Chain} with 4 chains."
    show(buff, MIME"text/plain"(), last(eachchain(pdb)))
    @test length(split(String(take!(buff)))) == 14*7 + 8

end