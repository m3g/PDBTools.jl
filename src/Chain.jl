
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
        error("Range $range does not correspond to a single residue or molecule.")
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
    if i <= 0 || i > length(chain)
        throw(ArgumentError("Index must be in 1:$(length(chain)). Attempted to access index $i."))
    end
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

@testitem "Chains" begin
    ats = read_pdb(PDBTools.TESTPDB)
    chains = eachchain(ats)
    @test length(chains) == 1
end


#atoms = readPDB("oligomer(1).pdb")
#length(atoms)
#atoms[1066].chain
#Chain(atoms, 1:1065)
#Residue(atoms, 1:13)

#rx = for atoms in first(eachresidue(atoms))
#    println(atoms)
#end
#xs = eachchain(atoms)
#length(xs)

#chainss =  collect(eachchain(atoms))

#x = for chains in eachchain(atoms)
#    println(chains)
#    end

#y = for res in eachresidue(atoms)
#    println(res)
#end