#
# Structure to read mmCIF field data indexes
#
Base.@kwdef mutable struct Indexes_mmCIF_fields
    index::Int = 0
    name::Int = 0
    resname::Int = 0
    chain::Int = 0
    resnum::Int = 0
    x::Int = 0
    y::Int = 0
    z::Int = 0
    beta::Int = 0
    occup::Int = 0
end
