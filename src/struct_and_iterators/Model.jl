"""
    Model

Model data structure. It carries the data of a model in a PDB file.
Models must be consecutive in the `atoms` vector, and
are identified by having the same `model` field.
 
The Model structure carries the properties of the model
it contains, but it does not copy the original vector of atoms, only the model
meta data and the reference to the original vector. Thus, changes
in the model atoms will be reflected in the original vector of atoms.

### Example

In the example below, 8S8N is PDB entry with 11 models.

```jldoctest
julia> using PDBTools

julia> ats = wget("8S8N");

julia> models = collect(eachmodel(ats))
11-element Vector{Model}[
    1-(234 atoms))
    2-(234 atoms))
    ⋮
    10-(234 atoms))
    11-(234 atoms))
]

julia> models[1]
 Model 1 with 234 atoms.
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     DLE     A        2        1   -5.811   -0.380   -2.159  1.00  0.00     1                 1
       2   CA     DLE     A        2        1   -4.785   -0.493   -3.227  1.00  0.00     1                 2
⋮
     233  HT2   A1H5T     B      101       13   -5.695    5.959   -3.901  1.00  0.00     1               233
     234  HT1   A1H5T     B      101       13   -4.693    4.974   -2.743  1.00  0.00     1               234

```

"""
struct Model{T<:Atom,Vec<:AbstractVector{T}} <: AbstractStructuralElement{T}
    atoms::Vec
    range::UnitRange{Int}
    number::Int32
end

# Necessary for the interface: define the same_struct_element function
same_struct_element(::Type{Model}, at1::Atom, at2::Atom) = at1.model == at2.model

# Constructors
function Model(atoms::AbstractVector{<:Atom}, range::AbstractRange{<:Integer})
    _check_unique(Model, atoms, range)
    return Model(atoms, UnitRange{Int}(range), model(first(atoms[range])))
end
Model(atoms::AbstractVector{<:Atom}) = Model(atoms, eachindex(atoms))

#
# Define the iterator
#
"""
    eachmodel(atoms::AbstractVector{<:Atom})

Iterator for the models of a selection.

### Example

Here we show how to iterate over the models of a PDB file, annotate 
the index of the first atom of each model, and collect all models.

```jldoctest
julia> using PDBTools

julia> ats = wget("8S8N");

julia> models = eachmodel(ats)
 Model iterator with length = 11

julia> first_atom = Atom[]
       for model in models
           push!(first_atom, model[1])
       end
       @show index.(first_atom);
index.(first_atom) = Int32[1, 235, 469, 703, 937, 1171, 1405, 1639, 1873, 2107, 2341]

julia> collect(models)
11-element Vector{Model}[
    1-(234 atoms))
    2-(234 atoms))
    ⋮
    10-(234 atoms))
    11-(234 atoms))
]
```

"""
eachmodel(atoms::AbstractVector{<:Atom}) = EachStructuralElement{Model}(atoms)

# Specific getters for this type
model(model::Model) = model.number

#
# Testing interface
#
@testitem "Model iterator" begin
    using PDBTools
    atoms = wget("8S8N")
    models = eachmodel(atoms)
    @test length(models) == 11
    @test firstindex(models) == 1
    @test lastindex(models) == 11
    @test length(last(models)) == 234
    @test_throws ArgumentError Model(atoms, 230:240)
    @test_throws ArgumentError models[1]
    m = collect(models)
    @test model(m[1]) == 1
    @test model(m[2]) == 2
    @test m[1].range == 1:234
    @test m[2].range == 235:468
    @test length(m[1]) == 234
    @test mass(m[1]) ≈ 1561.8828799999928
    @test size(m[1]) == (234,)
    @test eltype(m[1]) == Atom
    @test sum(mass(at) for at in m[1]) ≈ 1561.8828799999928
    m = Model(atoms[1:234])
    @test model(m) == 1
end

#
# io show functions
#
function Base.show(io::IO, mod::Model)
    compact = get(io, :compact, false)::Bool
    if compact
        print(io, "$(model(mod))-($(length(mod)) atoms))")
    else
        println(io, " Model $(model(mod)) with $(length(mod)) atoms.")
        show(IOContext(io, :type => false), @view mod.atoms[mod.range])
    end
end

@testitem "Model show" begin
    using PDBTools
    using ShowMethodTesting
    ENV["LINES"] = 10
    ENV["COLUMNS"] = 120
    ats = wget("8S8N")
    m = eachmodel(ats)
    @test parse_show(m; repl=Dict("PDBTools." => "")) ≈ "Model iterator with length = 11"
    mc = collect(m)
    @test parse_show(mc; repl=Dict("PDBTools." => "")) ≈ """
    11-element Vector{Model}[ 
    1-(234 atoms))
    2-(234 atoms))
    ⋮
    10-(234 atoms))
    11-(234 atoms))
    """
    @test parse_show(mc[1]; repl=Dict(r"^((?:[^\n]*\n){3}).*"s => s"\1")) ≈ """
     Model 1 with 234 atoms.
     index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     DLE     A        2        1   -5.811   -0.380   -2.159  1.00  0.00     1       1         1
       2   CA     DLE     A        2        1   -4.785   -0.493   -3.227  1.00  0.00     1       1         2
    """
end
