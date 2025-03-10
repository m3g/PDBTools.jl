abstract type AbstractStructuralElement{T} <:AbstractVector{T} end

mass(s::AbstractStructuralElement) = mass(@view s.atoms[s.range])

function Base.getindex(s::AbstractStructuralElement, i::Integer)
    i >= 0 || throw(ArgumentError("Index must be in 1:$(length(s))"))
    (i <= length(s)) || throw(ArgumentError("$(typeof(s)) has $(length(model)) atoms, tried to fetch index $i."))
    i = first(s.range) + i - 1
    s.atoms[i]
end

#
# Structure and function to define the iterators
#
struct EachStructuralElement{STYPE, T<:AbstractVector{<:Atom}}
    atoms::T
end
EachStructuralElement{STYPE}(atoms) where {STYPE} = EachStructuralElement{STYPE, typeof(atoms)}(atoms)

Base.collect(sit::EachStructuralElement{STYPE}) where {STYPE} = collect(STYPE, sit)
Base.length(sit::EachStructuralElement) = sum(1 for s in sit)
Base.firstindex(sit::EachStructuralElement) = 1
Base.lastindex(sit::EachStructuralElement) = length(sit)
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
_same(type::Type{AbstractStructuralElement}, ::Atom, ::Atom) = error("_same internal interface not implemented for $(type)")
function Base.iterate(sit::EachStructuralElement{STYPE}, current_atom=firstindex(sit.atoms)) where {STYPE}
    current_atom > length(sit.atoms) && return nothing
    next_atom = current_atom + 1
    while next_atom <= length(sit.atoms) && _same(STYPE, sit.atoms[current_atom], sit.atoms[next_atom])
        next_atom += 1
    end
    return (STYPE(sit.atoms, current_atom:next_atom-1), next_atom)
end

#
# Iterate over the structural element
#
function Base.iterate(s::AbstractStructuralElement, current_atom=firstindex(s))
    current_atom > lastindex(s) && return nothing
    return (s[current_atom], current_atom + 1)
end
