""" 

###Examples
```julia-repl
julia> pdb = read_pdb("protein_test.pdb")

julia> for chains in eachchain(pdb)
    println("  Chain: $(name(chains))")
    println(" with $(length(collect(eachresidue(chains)))) residues")
    println("each chain has $(length(chains)) atoms")
    end
  Chain: A
 with 3 residues
each chain has 48 atoms
  Chain: B
 with 3 residues
each chain has 48 atoms
  Chain: C
 with 3 residues
each chain has 48 atoms

for chain in eachchain(pdb)
    println(" Chain: $(name(chain))")
    println("$(length(eachresidue(chain))) residues") 

    for res in eachresidue(chain)
        println("$(name(res))")
        println(" with $(length(res)) atoms")
    end
end
 Chain: A
3 residues
ASP
 with 12 atoms
GLN
 with 17 atoms
LEU
 with 19 atoms
 Chain: B
3 residues
ASP
 with 12 atoms
GLN
 with 17 atoms
LEU
 with 19 atoms
 Chain: C
3 residues
ASP
 with 12 atoms
GLN
 with 17 atoms
LEU
 with 19 atoms

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

@testitem "Chain iterator" begin
    pdb = read_pdb("protein_test.pdb")
    chains = eachchain(pdb)
    @test Chain(pdb, 1:48).range == 1:48
    @test length(chains) == 3
    @test firstindex(chains) == 1
    @test lastindex(chains) == 3
    @test_throws ArgumentError chains[1]
    chains = collect(eachchain(pdb))
    @test name(chains[3]) == "C"
    @test index.(filter(at -> resname(at) == "ASP" && name(at) == "CA", chains[1])) == [2]
    @test length(findall(at -> resname(at) == "GLN", chains[1])) == 17
end

