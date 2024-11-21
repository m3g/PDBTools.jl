#
# Formula data type. Contains the number of atoms of each type in a vector of tuples.
#
struct Formula
    formula::Vector{Tuple{Atom,Int}}
end
Base.getindex(f::Formula, i) = (String(pdb_element(first(f.formula[i]))), last(f.formula[i]))

const sub_int = (
    "0" => "₀",
    "1" => "₁",
    "2" => "₂",
    "3" => "₃",
    "4" => "₄",
    "5" => "₅",
    "6" => "₆",
    "7" => "₇",
    "8" => "₈",
    "9" => "₉",
)

"""
    formula(atoms::AbstractVector{<:Atom})

Returns the molecular formula of the current selection. The output is an indexable
"Formula" structure, where each element is a tuple with the element name and the number of atoms.

## Example

```jldoctest
julia> using PDBTools

julia> pdb  = read_pdb(PDBTools.TESTPDB, "residue 1"); # testing PDB file

julia> resname(pdb[1])
"ALA"

julia> f = formula(pdb)
H₇C₃N₁O₁

julia> f[1]
("H", 7)

```

"""
function formula(atoms::AbstractVector{<:AtomType}) where {AtomType<:Atom}
    f = Formula(Tuple{AtomType,Int}[])
    for at in atoms
        i = findfirst(el -> pdb_element(first(el)) == element(at), f.formula)
        if isnothing(i)
            push!(f.formula, (Atom(pdb_element=element(at)), 1))
        else
            f.formula[i] = (first(f.formula[i]), last(f.formula[i]) + 1)
        end
    end
    # Sort by atomic number
    sort!(f.formula, by = el -> atomic_number(first(el)))
    return f
end
function Base.show(io::IO, f::Formula)
    s = ""
    for el in f.formula
        s *= pdb_element(first(el)) * format(last(el))
    end
    for sub in sub_int
        s = replace(s, sub)
    end
    print(io, s)
end

"""
    stoichiometry(atoms::AbstractVector{<:Atom})

Returns the stoichiometry of atom selection in a `Formula` structure. 

### Example

```jldoctest
julia> using PDBTools

julia> pdb  = read_pdb(PDBTools.TESTPDB, "water"); # testing PDB file

julia> stoichiometry(pdb)
H₂O₁
```

"""
function stoichiometry(atoms::AbstractVector{<:Atom})
    f = formula(atoms)
    d = gcd((x[2] for x in f.formula)...)
    for (i, p) in pairs(f.formula)
        f.formula[i] = (p[1], p[2] ÷ d)
    end
    f
end

@testitem "Formula" begin
    using PDBTools
    pdb  = read_pdb(PDBTools.TESTPDB, "resname GLY")
    f = formula(pdb)
    buff = IOBuffer()
    show(buff, f)
    @test String(take!(buff)) == "H₃₆C₂₄N₁₂O₁₂"
    @test f[1] == ("H", 36)
    @test f[4] == ("O", 12)
    s = stoichiometry(pdb)
    buff = IOBuffer()
    show(buff, s)
    @test String(take!(buff)) == "H₃C₂N₁O₁"
    @test s[1] == ("H", 3)
    @test s[4] == ("O", 1)
end
