
_tryparse(::Type{T}, s) where {T<:AbstractFloat} = tryparse(T, s)
function _tryparse(::Type{T}, s) where {T<:Integer} 
    i = tryparse(T, s) 
    # try to parse hexadecimal number
    if isnothing(i) 
        i = tryparse(T, s, base=16)
    end
    return i
end

function _parse(
    ::Type{T}, 
    string, 
    range=firstindex(string):lastindex(string); 
    alt::Union{T,Nothing}=nothing
) where {T<:Real}
    s = @view(string[firstindex(range):min(lastindex(range),lastindex(string))])
    x = _tryparse(T, s)
    !isnothing(x) && return x
    if isnothing(alt)
        throw(ArgumentError("Could not read $T from string: \"$s\""))
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
    isempty(range) && return alt
    first_char = if isspace(string[first(range)])
        findfirst(!isspace, string)
    else
        first(range)
    end
    last_char =  if isspace(string[last(range)])
        findlast(!isspace, string)
    else
        last(range)
    end
    if isnothing(first_char) | isnothing(last_char)
        return alt
    end
    return S(@view(string[first_char:last_char]))
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

function _parse(::Type{T}, string, range; alt=nothing) where {T<:Nothing}
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

## Alternative implementation using generated functions (same peformance as far as tested)
## https://discourse.julialang.org/t/unroll-setfield/122545/22?u=lmiq
#@generated function setfield_generated!(atom, field_values::FIELDS, inds_and_names::TUPTUP, ::Val{N}) where {FIELDS,TUPTUP,N}
#    quote
#        @inline
#        Base.@nexprs $N i -> begin
#            ifield, valfield = inds_and_names[i]
#            field = unwrap(valfield)
#            T = typeof(getfield(atom, field))
#            setfield!(atom, field, _parse(T, field_values[ifield]; alt=_alt(T)))
#        end
#    end
#end

@testitem "_parse" begin
    using PDBTools: _parse, _parse_charge
    @test _parse(Int, "  1  ") == 1
    @test _parse(Int, "  A  ") == 10
    @test _parse(Float32, "  1.0  ") == 1.0f0
    @test_throws ArgumentError _parse(Int, "  ???  ") 
    @test_throws ArgumentError _parse(Float32, "  A  ") 
    @test _parse(String, "  A  ") == "A"
    @test _parse(String, "") === nothing
    @test _parse_charge("1+") == "1"
    @test _parse_charge("+1") == "1"
    @test _parse_charge("1-") == "-1"
    @test _parse_charge("-1") == "-1"
end
