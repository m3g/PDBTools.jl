export ContactMap, contact_map, residue_residue_distance

"""
    residue_residue_distance(
        r1::PDBTools.Residue, 
        r2::PDBTools.Residue; 
        positions::AbstractVector{AbstractVector{T}}=nothing; 
        unitcell=nothing
    ) 

Calculate the minimum distance between two residues in a protein structure. 
If the `positions` argument is not provided, the function calculates the distance
using the coordinates of the atoms in the residues. If `positions` is provided,
the function uses the coordinates in the positions array. 

# Arguments

- `r1::PDBTools.Residue`: Residue 1
- `r2::PDBTools.Residue`: Residue 2
- `positions::AbstractVector{AbstractVector{T}}`: Optional alternate positions of the atoms in the structure.
- `unitcell=nothing`: Optional unit cell dimensions for periodic boundary conditions.

!!! note
    The index of the atoms in the residues must match the index of the atoms in the
    positions array. 

# Example

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> residues = collect(eachresidue(ats));

julia> r1 = residues[1]; r10 = residues[10];

julia> println(name(r1), resnum(r1), " and ", name(r10), resnum(r10))
LYS211 and GLU220

julia> d = residue_residue_distance(r1, r10)
16.16511f0
```

"""
function residue_residue_distance(
    r1::PDBTools.Residue,
    r2::PDBTools.Residue;
    positions::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing,
    unitcell=nothing
)
    dmin = typemax(first(r1).x)
    for (iat, jat) in Iterators.product(eachindex(r1), eachindex(r2))
        p_i = isnothing(positions) ? PDBTools.coor(r1[iat]) : positions[PDBTools.index(r1[iat])]
        p_j = isnothing(positions) ? PDBTools.coor(r2[jat]) : positions[PDBTools.index(r2[jat])]
        p_j = !isnothing(unitcell) ? wrap(p_j, p_i, unitcell) : p_j
        d = norm(p_j - p_i)
        d < dmin && (dmin = d)
    end
    return dmin
end

@testitem "residue_residue_distance" begin
    using PDBTools
    ats = read_pdb(PDBTools.TESTPDB, "protein")
    residues = collect(eachresidue(ats))
    r1 = residues[1]
    r10 = residues[10]
    # Testing call with residue information only
    @test residue_residue_distance(r1, r10) ≈ 5.6703672f0
    d = residue_residue_distance(r1, r10; positions=coor(ats))
    @test d ≈ 5.6703672f0

    # Test with a PBC cell (protein broken by the PBCs)
    pdb_pbc = read_pdb(PDBTools.TESTPBC, "protein")
    uc = read_unitcell(PDBTools.TESTPBC)
    pdb_nopbc = read_pdb(PDBTools.TESTNOPBC, "protein")
    r_pbc = collect(eachresidue(pdb_pbc))
    r_nopbc = collect(eachresidue(pdb_nopbc))
    for i in eachindex(r_pbc, r_nopbc)
        for j in eachindex(r_pbc, r_nopbc)
            dpbc = residue_residue_distance(r_pbc[i], r_pbc[j]; unitcell=uc)
            dnopbc = residue_residue_distance(r_nopbc[i], r_nopbc[j])
            @test dpbc ≈ dnopbc rtol=1e-2
        end
    end
end

"""
    ContactMap{Bool|Real}

Data structure to store contact maps between residues in a protein structure.
The contact map is a matrix of distances between residues. A contact is defined 
when the distance between any two atoms of the residues is less than a given threshold `dmax`.

If the distance between two residues is greater than `dmax`, the value in the
matrix is `missing`, indicating that there is no contact between the residues.
If the distance is less than `dmax`, the value in the matrix is the distance
between the residues.

The `gap` parameter is used to calculate contacts between residues separated by
a given number of residues. For example, if `gap=3`, the contact map was 
calculated between residues separated by at least 3 residues in the sequence.

# Fields

- `matrix::Matrix{Union{Missing,T}}`: Matrix of distances between residues.
- `d::T`: Threshold distance for a contact.
- `gap::Int`: Gap between residues to calculate contacts.

If the contact map was calculated with `discrete=true`, the matrix contains
`Bool` values, where `true` indicates a contact and `false` indicates no contact.
On the other hand, if `discrete=false`, the matrix contains distances between
residues.

"""
struct ContactMap{T<:Union{Union{Missing,Bool},Union{Missing,<:Real}}}
    matrix::Matrix{T}
    d::Float64
    gap::Int
    residues1::Vector{PDBTools.Residue}
    residues2::Vector{PDBTools.Residue}
end
function Base.zero(
    ::Type{ContactMap{T}},
    n, m, d::Real, gap::Int,
    residues1::Vector{PDBTools.Residue},
    residues2::Vector{PDBTools.Residue},
) where {T}
    ContactMap{T}([missing for i in 1:n, j in 1:m], d, gap, residues1, residues2)
end
Base.setindex!(map::ContactMap, value, i, j) = map.matrix[i, j] = value
Base.getindex(map::ContactMap, i, j) = map.matrix[i, j]
function Base.show(io::IO, ::MIME"text/plain", map::ContactMap{T}) where {T}
    print(io, "ContactMap{$T} of size $(size(map.matrix)), with threshold $(map.d) and gap $(map.gap)")
end

"""
    contact_map(
        atoms1::AbstractVector{<:PDBTools.Atom}
        atoms2::AbstractVector{<:PDBTools.Atom}; # optional for contacts between two structures
        dmax::Real=4.0,
        gap::Int=0, # only available if atoms2 is not provided
        unitcell=nothing,
        discrete::Bool=true,
        positions::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing,
    )

Calculate the contact map between residues in a protein* structure (if only 
`atoms1` is provided) or between residues in two different protein* structures
(`atoms1` and `atoms2`). 

!!! note 
    The distance used to define a contact is the **minimum** distance between
    any two atoms of the residues of the atoms groups provided, with a 
    threshold distance `dmax`.

Returns the contact map as a `ContactMap` object.

*The term "protein" is used here to refer to any structure with residues.

# Arguments

- `atoms1::AbstractVector{<:PDBTools.Atom}`: Atoms of the first structure.
- `atoms2::AbstractVector{<:PDBTools.Atom}`: Atoms of the second structure, if provided. 

# Optional keyword arguments

- `dmax::Real=4.0`: Threshold distance for a contact.
- `gap::Int=0`: Gap between residues to calculate contacts.
- `unitcell=nothing`: Unit cell dimensions for periodic boundary conditions.
- `discrete::Bool=true`: If `true`, the matrix contains `Bool` values, where `true` 
  indicates a contact and `false` indicates no contact. 
  If `false`, the matrix contains distances between residues.
- `positions`: Positions of the atoms in the structure. If provided, the function uses these 
  positions to calculate the distance between residues.

# Examples

## Contact map between residues in the same structure

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> cA = select(ats, "chain A");

julia> cB = select(ats, "chain B");

julia> map = contact_map(cA, cB) # contact map between chains A and B
ContactMap{Union{Missing, Bool}} of size (243, 12), with threshold 4.0 and gap 0 

julia> # using Plots; heatmap(map); # uncomment to plot the contact map
```

## Contact map between residues in two different structures

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> cA = select(ats, "chain A");

julia> cB = select(ats, "chain B");

julia> map = contact_map(cA, cB) # contact map between chains A and B
ContactMap{Union{Missing, Bool}} of size (243, 12), with threshold 4.0 and gap 0 

julia> # using Plots; heatmap(map); # uncomment plot the contact map
```

"""
function contact_map end

function contact_map(
    atoms1::AbstractVector{<:PDBTools.Atom};
    dmax::Real=4.0,
    gap::Int=0,
    discrete::Bool=true,
    unitcell=nothing,
    positions::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing,
)
    type = Union{Missing,discrete ? Bool : typeof(atoms1[1].x)}
    residues = collect(PDBTools.eachresidue(atoms1))
    map = zero(ContactMap{type}, length(residues), length(residues), dmax, gap, residues, residues)
    for ires in eachindex(residues), jres in ires+gap:length(residues)
        r1 = residues[ires]
        r2 = residues[jres]
        d12 = residue_residue_distance(r1, r2; positions, unitcell)
        if discrete
            map[ires, jres] = d12 <= dmax ? true : false
        else
            map[ires, jres] = d12 <= dmax ? d12 : missing
        end
        map[jres, ires] = map.matrix[ires, jres]
    end
    return map
end

function contact_map(
    atoms1::AbstractVector{<:PDBTools.Atom},
    atoms2::AbstractVector{<:PDBTools.Atom};
    dmax::Real=4.0,
    discrete::Bool=true,
    unitcell=nothing,
    positions::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing,
)
    type = Union{Missing,discrete ? Bool : typeof(atoms1[1].x)}
    residues = collect(PDBTools.eachresidue(atoms1))
    residues2 = collect(PDBTools.eachresidue(atoms2))
    map = zero(ContactMap{type}, length(residues), length(residues2), dmax, 0, residues, residues2)
    for ires in eachindex(residues), jres in eachindex(residues2)
        r1 = residues[ires]
        r2 = residues2[jres]
        d12 = residue_residue_distance(r1, r2; positions, unitcell)
        if discrete
            map[ires, jres] = d12 <= dmax ? true : false
        else
            map[ires, jres] = d12 <= dmax ? d12 : missing
        end
    end
    return map
end

@testitem "contact_map" begin
    using PDBTools
    # monomer
    ats = read_pdb(PDBTools.TESTPDB, "protein")
    map = contact_map(ats)
    @test size(map.matrix) == (104, 104)
    @test count(map.matrix) == 1106
    map = contact_map(ats; discrete=false)
    @test sum(skipmissing(map.matrix)) ≈ 2407.2163f0
    # dimer
    ats = read_pdb(PDBTools.DIMERPDB)
    cA = select(ats, "chain A")
    cB = select(ats, "chain B")
    map = contact_map(cA, cB)
    @test sum(map.matrix) == 17
    @test map[235, :] == [false, false, true, false, false, false, false, false, false, false, false, false]
    @test count(map[:, 3]) == 3
    map = contact_map(cA, cB; discrete=false)
    @test sum(skipmissing(map.matrix)) ≈ 58.00371f0

    # Test with and without periodic boundary conditions
    pdb_pbc = read_pdb(PDBTools.TESTPBC, "protein")
    uc = read_unitcell(PDBTools.TESTPBC)
    pdb_nopbc = read_pdb(PDBTools.TESTNOPBC, "protein")
    pbc = contact_map(pdb_pbc; discrete=true, unitcell=uc)
    no_pbc = contact_map(pdb_nopbc; discrete=true)
    @test pbc.matrix == no_pbc.matrix
    pbc = contact_map(pdb_pbc; discrete=false, unitcell=uc)
    no_pbc = contact_map(pdb_nopbc; discrete=false)
    @test all(ismissing(pbc.matrix[i]) | isapprox(pbc.matrix[i], no_pbc.matrix[i]; rtol=1e-2) for i in eachindex(pbc.matrix))

end