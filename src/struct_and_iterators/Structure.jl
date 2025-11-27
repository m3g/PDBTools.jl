"""
    Structure

A `Structure` is the main data structure in PDBTools.jl representing a molecular structure.
It is a wrapper around a vector of `Atom` objects, with additional metadata stored in a dictionary.

"""
struct Structure{V,T} <: AbstractVector{T}
    atoms::V
    _meta_data::Dict{Symbol,Any}
end
function Structure(ats::AbstractVector{<:Atom}=Atom{Nothing}[], _meta_data=Dict{Symbol,Any}(); kargs...) 
    s = Structure{typeof(ats), eltype(ats)}(ats, _meta_data)
    for (key, val) in pairs(kargs)
        s._meta_data[key] = val
    end
    return s
end
Base.length(s::Structure) = length(getfield(s, :atoms))
Base.size(s::Structure) = size(getfield(s, :atoms))
Base.getindex(s::Structure, i::Integer) = getfield(s, :atoms)[i]
Base.setindex!(s::Structure, at::Atom, i::Integer) = setindex!(getfield(s, :atoms), at, i)
Base.getindex(s::Structure, range::AbstractUnitRange{<:Integer}) = Structure(@view(s.atoms[range]), s._meta_data)

Base.setproperty!(s::Structure, field::Symbol, val) = setproperty!(s::Structure, Val(field), val)
Base.setproperty!(s::Structure, ::Val{field}, val) where {field} = getfield(s, :_meta_data)[field] = val
Base.setproperty!(s::Structure, ::Val{:atoms}, val::AbstractVector{<:Atom}) = setfield!(s, :atoms, val)
Base.setproperty!(_::Structure, ::Val{:_meta_data}, val) = throw(ArgumentError(" Cannot modify internal _meta_data field."))

Base.getproperty(s::Structure, ::Val{:atoms}) = getfield(s, :atoms)
Base.getproperty(s::Structure, ::Val{:_meta_data}) = getfield(s, :_meta_data)
Base.getproperty(s::Structure, field::Symbol) = getproperty(s, Val(field))
function Base.getproperty(s::Structure, ::Val{field}) where {field} 
    _meta_data = getfield(s, :_meta_data)
    return _meta_data[field]
end
Base.propertynames(s::Structure) = (:atoms, keys(s._meta_data)...)

function Base.show(io::IO, s::Structure)
    lines, cols = _displaysize(io)
    compact = get(io, :compact, false)::Bool
    ntitlelines = 6
    if !compact && cols >= 115 && lines > ntitlelines - 1 
        println(io, "   $(typeof(s))")
        println(io, "   number of atoms: $(length(s)), data: ", join(propertynames(s), ", "))
    end
    show(IOContext(io, :type => false, :ntitlelines => ntitlelines), s.atoms)
end

@testitem "Structure" begin
    using PDBTools
    s = Structure(read_pdb(PDBTools.TESTPDB))
    ats = s.atoms
    @test s[1] == ats[1]
    @test s[end] == ats[end]
    @test s[5:10] == ats[5:10]
    @test size(s) == size(ats)
    @test length(s) == length(ats)
    s[2] = Atom(name="CA")
    @test s[2] == Atom(name="CA")
    @test propertynames(s) == (:atoms,)
    s.unitcell = [ 1 0 0 ; 0 1 0 ; 0 0 1]
    @test s.unitcell == [ 1 0 0 ; 0 1 0 ; 0 0 1]
    @test propertynames(s) == (:atoms, :unitcell)
    s = Structure(ats; unitcell=[1 0 0; 0 1 0; 0 0 1])
    @test propertynames(s) == (:atoms, :unitcell)
    @test s.unitcell == [ 1 0 0 ; 0 1 0 ; 0 0 1]
    @test_throws "Cannot modify" s._meta_data = 1
end

@testitem "Structure - show" begin
    using PDBTools
    using ShowMethodTesting
    ENV["LINES"] = 10
    ENV["COLUMNS"] = 120
    at = Atom(name="CA", segname="X")
    s = Structure([at]; unitcell=[1 0 0; 0 1 0; 0 0 1])
    @test parse_show(s; repl=Dict("PDBTools." => "")) ≈ """
       Structure{Vector{Atom{Nothing}}, Atom{Nothing}}
       number of atoms: 1, data: atoms, unitcell
       index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       0   CA     XXX     X        0        0    0.000    0.000    0.000  0.00  0.00     0       X         0
    """
    @test parse_show(Structure([at for _ in 1:50]); repl=Dict(r"^((?:[^\n]*\n){4}).*"s => s"\1", r"PDBTools." => "")) ≈ """
       Structure{Vector{Atom{Nothing}}, Atom{Nothing}}
   number of atoms: 50, data: atoms
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       0   CA     XXX     X        0        0    0.000    0.000    0.000  0.00  0.00     0       X         0
    """
    @test parse_show(Dict(1 => Structure([at, at, at])); repl=Dict("PDBTools." => "")) ≈ """
    Dict{Int64, Structure{Vector{Atom{Nothing}}, Atom{Nothing}}} with 1 entry:
    1 => [ Atom(0CA-XXX0X), Atom(0CA-XXX0X), Atom(0CA-XXX0X) ]
    """
    @test parse_show([Structure([at, at]), Structure([at, at, at])]; repl=Dict("PDBTools." => "")) ≈ """
    2-element Vector{Structure{Vector{Atom{Nothing}}, Atom{Nothing}}}[ 
    [ Atom(0CA-XXX0X), Atom(0CA-XXX0X) ]
    [ Atom(0CA-XXX0X), Atom(0CA-XXX0X), Atom(0CA-XXX0X) ]
    ]
    """
    @test parse_show(Atom{Nothing}[]; repl=Dict("PDBTools." => "")) ≈ """
       Vector{Atom{Nothing}} with 0 atoms with fields:
    index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
    """
    @test parse_show(Structure(); repl=Dict("PDBTools." => "")) ≈ """
       Structure{Vector{Atom{Nothing}}, Atom{Nothing}}
    number of atoms: 0, data: atoms
    index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
    """
    @test parse_show([ Structure() ]; repl=Dict("PDBTools." => "")) ≈ """
    1-element Vector{Structure{Vector{Atom{Nothing}}, Atom{Nothing}}}[ 
        [  ]
    ]
    """
    @test parse_show(Structure([ at for _ in 1:20 ]); mime=nothing, context = :compact => true) ≈ """
    [ Atom( 0 CA-XXX 0 X), Atom( 0 CA-XXX 0 X), Atom( 0 CA-XXX 0 X), Atom( 0 CA-XXX 0 X)…
    """
end