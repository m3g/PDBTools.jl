function _parse(::Type{T}, string, range=nothing) where {T<:AbstractFloat}
    s = if isnothing(range)
        string
    else
        @view(string[range[begin]:min(range[end],length(string))])
    end
    x = tryparse(T, s)
    return isnothing(x) ? zero(T) : x
end

function _parse(::Type{T}, string, range=nothing; alt::Union{T,Nothing}=nothing) where {T<:Integer}
    s = if isnothing(range)
        string
    else
        @view(string[range[begin]:min(range[end],length(string))])
    end
    i = tryparse(T, s)
    !isnothing(i) && return i
    # try to parse as hexadecimal number
    i = tryparse(T, s, base=16)
    !isnothing(i) && return i
    if isnothing(alt)
        error("Could not read integer from string: \"$s\"")
    else
        return alt
    end
end

function _parse(::Type{S}, string, range=nothing; alt::Union{S,Nothing}=nothing) where {S<:AbstractString}
    s, range = if isnothing(range)
        string, 1:length(string)
    else
        @view(string[range[begin]:min(range[end],length(string))]), range
    end
    length(range) > 0 || return isnothing(alt) ? S(" ") : alt
    first_char = findfirst(>(' '), s)
    isnothing(first_char) && return isnothing(alt) ? S(" ") : alt
    first_char = first(range) + first_char - 1
    last_char = first(range) + findlast(>(' '), s) - 1
    return S(string[first_char:last_char])
end

function _parse(::Type{T}, string, range=1:1) where {T<:Nothing}
    return nothing
end

#
# Great contribution from Mason Protter
# https://discourse.julialang.org/t/unroll-setfield/122545/3?u=lmiq 
#
#=

field_values is a tuple of values to set the fields of atom to.
inds_and_names is a tuple of tuples of indices and field names to set the fields of atom to.

ex: 

field_values = (1, 2.0, "CA")
inds_and_names = ((1, Val(:index)), (2, Val(:occup)), (3, Val(:name)))
setfield_recursive!(atom, field_values, inds_and_names)

would set atom.index = 1, atom.occup = 2.0, atom.name = "CA"

=#
unwrap(::Val{T}) where {T} = T
function setfield_recursive!(atom::AtomType, field_values::FIELDS, inds_and_names::TUPTUP) where {AtomType, FIELDS, TUPTUP}
    isempty(inds_and_names) && return atom
    i, valfield = first(inds_and_names)
    field = unwrap(valfield)
    T = typeof(getfield(atom, field))
    setfield!(atom, field, _parse(T, field_values[i]))
    setfield_recursive!(atom, field_values, Base.tail(inds_and_names))
end
