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
3-element Vector{InlineStrings.String3}:
 "A"
 "B"
 "C"

julia> chains = collect(eachchain(ats))
   Array{Chain,1} with 3 chains.

julia> chains[2].range
49:96

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
    i > 0 || throw(ArgumentError("Chain index must be in 1:$(length(residue))"))
    (i <= length(residue)) || throw(ArgumentError("Chain has $(length(residue)) atoms, tried to fetch index $i."))
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
3

julia> for chain in eachchain(ats)
              println(chain.range)
              end
1:48
49:96
97:144

```

"""
eachchain(atoms::AbstractVector{<:Atom}) = EachChain(atoms)

# Collect chains default constructor
Base.collect(c::EachChain) = collect(Chain, c)
Base.length(chains::EachChain) = sum(1 for chain in chains)
Base.firstindex(chains::EachChain) = 1
Base.lastindex(chains::EachChain) = length(chains)

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
    print(io, "   Array{Chain,1} with $(length(chains)) chains.")
end

@testitem "Chain iterator" begin
    pdb = read_pdb(PDBTools.CHAINSPDB)
    chains = eachchain(pdb)
    @test Chain(pdb, 1:48).range == 1:48
    @test_throws ArgumentError Chain(pdb, 49:97).range 
    @test length(chains) == 3
    @test firstindex(chains) == 1
    @test lastindex(chains) == 3
    @test_throws ArgumentError chains[1]
    chains = collect(eachchain(pdb))
    @test name(chains[3]) == "C"
    @test index.(filter(at -> resname(at) == "ASP" && name(at) == "CA", chains[1])) == [2]
    @test length(findall(at -> resname(at) == "GLN", chains[1])) == 17
    @test mass(chains[1]) == 353.37881000000016
    @test chains[3].segname == "ASYN"
    @test chains[2].model == 1
    @test chains[1].chain == "A"
    
end

