#
# Formula data type. Contains the number of atoms of each type in a vector of tuples.
#
struct Formula
    formula::Vector{Tuple{Atom,Int}}
end
Base.getindex(f::Formula, i) = f.formula[i]

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

Returns the molecular formula of the current selection. 

## Example

```jldoctest
julia> using PDBTools

julia> pdb  = read_pdb(PDBTools.TESTPDB, "residue 1"); # testing PDB file

julia> resname(pdb[1])
"ALA"

julia> formula(pdb)
H₇C₃N₁O₁
```

"""
function formula(atoms::AbstractVector{<:Atom})
    f = Formula(Tuple{Atom,Int}[])
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

```julia-repl
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
