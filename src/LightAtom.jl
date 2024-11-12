@kwdef mutable struct LightAtom <: AbstractAtom
    index::Int32 = 0
    index_pdb::Int32 = 0
    name::String7 = "X"
    resname::String7 = "XXX"
    chain::String3 = "X"
    resnum::Int32 = 0
    residue::Int32 = 0
    x::Float32 = 0.0
    y::Float32 = 0.0
    z::Float32 = 0.0
    beta::Float32 = 0.0
    occup::Float32 = 0.0
    charge::Float32 = 0.0
    model::Int32 = 0
    segname::String7 = ""
end