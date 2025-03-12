abstract type AbstractStructuralElement{T} <: AbstractVector{T} end

function Base.getindex(s::AbstractStructuralElement, i::Integer)
    (1 <= i <= length(s)) || throw(ArgumentError("$(typeof(s)) has $(length(s)) atoms, tried to fetch index $i."))
    i = first(s.range) + i - 1
    s.atoms[i]
end

# Necessary for the interface: define the same_struct_element function
same_struct_element(type::Type{<:AbstractStructuralElement}, ::Atom, ::Atom) =
    throw(ArgumentError("same_struct_element internal interface not implemented for $(type)"))

# Function that performs checks on the range
function _check_unique(::Type{T}, atoms, range) where {T<:AbstractStructuralElement}
    at1 = first(atoms[range])
    if any(!same_struct_element(T, at1, at2) for at2 in @view(atoms[range]))
        throw(ArgumentError("""\n 
                Range $range does not correspond to a single $T structural elements. 
                
        """))
    end
end

#
# Structure and function to define the iterators
#
struct EachStructuralElement{STYPE,T<:AbstractVector{<:Atom}}
    atoms::T
end
EachStructuralElement{STYPE}(atoms) where {STYPE} = EachStructuralElement{STYPE,typeof(atoms)}(atoms)

Base.collect(sit::EachStructuralElement{STYPE}) where {STYPE} = collect(STYPE, sit)
Base.length(sit::EachStructuralElement) = count(el -> true, sit)
Base.firstindex(sit::EachStructuralElement) = 1
Base.lastindex(sit::EachStructuralElement) = length(sit)
function Base.last(sit::EachStructuralElement{STYPE}) where {STYPE}
    atoms = sit.atoms
    ilast = lastindex(atoms)
    iprev = ilast - 1
    while iprev >= firstindex(atoms) && same_struct_element(STYPE, atoms[ilast], atoms[iprev])
        iprev -= 1
    end
    ifirst = iprev + 1
    return STYPE(atoms, ifirst:ilast)
end

function Base.getindex(::EachStructuralElement, ::Int)
    throw(ArgumentError("""\n
        The iterator does not support indexing. 
        Use collect(iterator) to get an indexable list.
    """))
end

# Array interface
Base.size(s::AbstractStructuralElement) = (length(s.range),)
Base.length(s::AbstractStructuralElement) = length(s.range)
Base.eltype(::AbstractStructuralElement) = Atom

#
# Iterate, lazily, over the iterator
#
function Base.iterate(sit::EachStructuralElement{STYPE}, current_atom=firstindex(sit.atoms)) where {STYPE}
    current_atom > length(sit.atoms) && return nothing
    next_atom = current_atom + 1
    while next_atom <= length(sit.atoms) && same_struct_element(STYPE, sit.atoms[current_atom], sit.atoms[next_atom])
        next_atom += 1
    end
    return (STYPE(sit.atoms, current_atom:next_atom-1), next_atom)
end

# Properties
mass(s::AbstractStructuralElement) = mass(@view s.atoms[s.range])

#
# Iterate over the structural element
#
function Base.iterate(s::AbstractStructuralElement, current_atom=firstindex(s))
    current_atom > lastindex(s) && return nothing
    return (s[current_atom], current_atom + 1)
end

#
# Default iterator show method
#
function Base.show(io::IO, sit::EachStructuralElement{STYPE}) where {STYPE}
    print(io, " $STYPE iterator with length = $(length(sit))")
end

@testitem "AbstractStructuralElement interface" begin
    using PDBTools
    struct Incomplete{T} <: PDBTools.AbstractStructuralElement{T} 
        atoms::Vector{T}
        range::UnitRange{Int}
    end
    function Incomplete(atoms::AbstractVector{<:Atom}, range::AbstractRange{<:Integer})
        PDBTools._check_unique(Incomplete, atoms, range)
    end
    Incomplete(atoms::AbstractVector{<:Atom}) = Incomplete(atoms, eachindex(atoms))
    ats = read_pdb(PDBTools.SMALLPDB)
    @test_throws ArgumentError Incomplete(ats)
end
