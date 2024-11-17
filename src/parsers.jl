function _parse(
    ::Type{T}, 
    string, 
    range=firstindex(string):lastindex(string); 
    alt=zero(T)
) where {T<:AbstractFloat}
    s = @view(string[firstindex(range):min(lastindex(range),lastindex(string))])
    x = tryparse(T, s)
    !isnothing(x) && return x
    if isnothing(alt)
        error("Could not read float from string: \"$s\"")
    else
        return alt
    end
end

function _parse(
    ::Type{T}, 
    string, 
    range=firstindex(string):lastindex(string); 
    alt::Union{T,Nothing}=nothing
) where {T<:Integer}
    s = @view(string[firstindex(range):min(lastindex(range),lastindex(string))])
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

function _parse(
    ::Type{S}, 
    string, 
    range=firstindex(string):lastindex(string); 
    alt::Union{S,Nothing}=nothing
) where {S<:AbstractString}
    s = @view(string[firstindex(range):min(lastindex(range),lastindex(string))])
    length(range) > 0 || return isnothing(alt) ? S("X") : alt
    first_char = findfirst(!isspace, s)
    isnothing(first_char) && return isnothing(alt) ? S("X") : alt
    last_char = findlast(!isspace, s)
    return S(s[first_char:last_char])
end

#
# In PDB, the charge sometimes appears as "+1" or "-1" and sometimes as "1+" or "1-"
#
function _parse_charge(charge::AbstractString)
    charge = if occursin('+', charge) 
        replace(charge, "+" => "")
    elseif occursin('-', charge)
        "-$(replace(charge, "-" => ""))"
    else
        charge
    end
    return charge
end

function _parse(::Type{T}, string, range=1:1; alt=nothing) where {T<:Nothing}
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
function _fast_setfield!(atom::AtomType, field_values::FIELDS, inds_and_names::TUPTUP) where {AtomType, FIELDS, TUPTUP} 
    setfield_recursive!(atom, field_values, inds_and_names)
# Alternative with generated function:
#    N = length(inds_and_names)
#    setfield_generated!(atom, field_values, inds_and_names, Val(N))
end

# Alternate values for fields that might be empty
_alt(::Type{S}) where {S<:AbstractString} = S("X")
_alt(::Type{T}) where {T} = zero(T)
# Unwrap Val-wrapped values
unwrap(::Val{T}) where {T} = T
function setfield_recursive!(atom::AtomType, field_values::FIELDS, inds_and_names::TUPTUP) where {AtomType, FIELDS, TUPTUP}
    isempty(inds_and_names) && return atom
    i, valfield = first(inds_and_names)
    field = unwrap(valfield)
    T = typeof(getfield(atom, field))
    setfield!(atom, field, _parse(T, field_values[i]; alt=_alt(T)))
    setfield_recursive!(atom, field_values, Base.tail(inds_and_names))
end

# Alternative implementation using generated functions (same peformance as far as tested)
# https://discourse.julialang.org/t/unroll-setfield/122545/22?u=lmiq
@generated function setfield_generated!(atom, field_values::FIELDS, inds_and_names::TUPTUP, ::Val{N}) where {FIELDS,TUPTUP,N}
    quote
        @inline
        Base.@nexprs $N i -> begin
            ifield, valfield = inds_and_names[i]
            field = unwrap(valfield)
            T = typeof(getfield(atom, field))
            setfield!(atom, field, _parse(T, field_values[ifield]; alt=_alt(T)))
        end
    end
end
