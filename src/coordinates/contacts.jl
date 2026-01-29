import CellListMap

export ContactMap, contact_map, residue_residue_distance

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

import Base: +, -

function _err_sum_cmap(c1::ContactMap, c2::ContactMap)
    if (length(c1.residues1) != length(c2.residues1)) ||
       (length(c1.residues2) != length(c2.residues2))
        throw(ArgumentError("""\n
            Arithmetic operations are only defined for contact maps with the same number of residues.

        """))
    end
    return nothing
end

function +(c1::ContactMap{Union{Missing,Bool}}, c2::ContactMap{Union{Missing,Bool}})
    _err_sum_cmap(c1, c2)
    c3_matrix = copy(c1.matrix)
    for i in eachindex(c1.matrix, c2.matrix, c3_matrix)
        c3_matrix[i] = if (c2.matrix[i] === true)
            true
        else
            c3_matrix[i]
        end
    end
    return ContactMap(c3_matrix, c1.d, c1.gap, c1.residues1, c1.residues2)
end
function -(c1::ContactMap{Union{Missing,Bool}}, c2::ContactMap{Union{Missing,Bool}})
    _err_sum_cmap(c1, c2)
    c3_matrix = zeros(Union{Int,Missing}, size(c1.matrix))
    for i in eachindex(c1.matrix, c2.matrix, c3_matrix)
        c3_matrix[i] = if (c1.matrix[i] === true) & (c2.matrix[i] === true)
            0
        elseif (c1.matrix[i] === true) & (!(c2.matrix[i] === true))
            1
        elseif (!(c1.matrix[i] === true)) & (c2.matrix[i] === true)
            -1
        else
            missing
        end
    end
    return ContactMap(c3_matrix, c1.d, c1.gap, c1.residues1, c1.residues2)
end

function +(c1::ContactMap{<:Union{Missing,<:AbstractFloat}}, c2::ContactMap{<:Union{Missing,<:AbstractFloat}})
    _err_sum_cmap(c1, c2)
    return ContactMap(c1.matrix + c2.matrix, c1.d, c1.gap, c1.residues1, c1.residues2)
end
function -(c1::ContactMap{<:Union{Missing,<:AbstractFloat}}, c2::ContactMap{<:Union{Missing,<:AbstractFloat}})
    _err_sum_cmap(c1, c2)
    return ContactMap(c2.matrix - c1.matrix, c1.d, c1.gap, c1.residues1, c1.residues2)
end

_err_same_type() = throw(ArgumentError("""\n
        Arithmetic operations are only defined for contact maps of the same type (discrete *or* continuous).
    """))
+(_::ContactMap, _::ContactMap) = _err_same_type()
-(_::ContactMap, _::ContactMap) = _err_same_type()

"""
    contact_map(
        atoms1::AbstractVector{<:PDBTools.Atom}
        atoms2::AbstractVector{<:PDBTools.Atom}; # optional for contacts between two structures
        dmax::Real=4.0,
        gap::Int=0, # only available if atoms2 is not provided
        unitcell=nothing,
        discrete::Bool=true,
        positions::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing,
        parallel::Bool=false,
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
- `parallel`: set to `true` to allow parallel computations. This is `false` by default because
  the calculation is usually fast and parallelization can use too much memory for very large systems.

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

#
# Fast contact map calculation for big structures using cell lists.
#
struct MapMatrix{T}
    m::Matrix{T}
end

_contact(x::Bool, y::Bool) = x | y
function _contact(x, y)
    ismissing(x) && ismissing(y) && return missing
    ismissing(x) && return y
    ismissing(y) && return x
    return min(x, y)
end

CellListMap.copy_output(x::MapMatrix) = MapMatrix(copy(x.m))
function CellListMap.reducer(x::MapMatrix, y::MapMatrix)
    for i in eachindex(x.m, y.m)
        x.m[i] = _contact(x.m[i], y.m[i])
    end
    return x
end
_reset(::Type{Bool}) = false
_reset(::Type{<:Real}) = missing
function CellListMap.reset_output!(x::MapMatrix{<:Union{Missing,<:T}}) where {T}
    x.m .= _reset(T)
    return x
end

function update_map_matrix(
    i, j, d2, map_matrix,
    index_residue::Vector{Int},
    gap::Real,
    discrete::Bool
)
    (; m) = map_matrix
    ires = index_residue[i]
    jres = index_residue[j]
    d_gap = abs(ires - jres)
    if d_gap >= gap
        if discrete
            m[ires, jres] = true
        elseif d_gap > 0 # same-residue distance is zero
            m[ires, jres] = _contact(m[ires, jres], sqrt(d2))
        end
        m[jres, ires] = m[ires, jres]
    end
    return map_matrix
end

function assign_index_residue(atoms, residues)
    residue_index = zeros(Int, length(atoms))
    iat = 0
    for (ires, res) in enumerate(residues)
        for at in res
            iat += 1
            residue_index[iat] = ires
        end
    end
    return residue_index
end

function contact_map(
    atoms1::AbstractVector{<:PDBTools.Atom};
    dmax::Real=4.0,
    gap::Int=0,
    discrete::Bool=true,
    unitcell=nothing,
    positions::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing,
    parallel::Bool=false,
)
    type = Union{Missing,discrete ? Bool : typeof(atoms1[1].x)}
    residues = collect(PDBTools.eachresidue(atoms1))
    index_residue = assign_index_residue(atoms1, residues)
    map = zero(ContactMap{type}, length(residues), length(residues), dmax, gap, residues, residues)
    map_matrix = MapMatrix(map.matrix)
    sys = CellListMap.ParticleSystem(
        positions=isnothing(positions) ? coor.(atoms1) : positions,
        cutoff=dmax,
        unitcell=unitcell,
        output=map_matrix,
        output_name=:map_matrix,
        parallel=parallel,
    )
    CellListMap.map_pairwise!(
        (x, y, i, j, d2, map_matrix) -> update_map_matrix(i, j, d2, map_matrix, index_residue, gap, discrete),
        sys
    )
    return map
end

function update_map_matrix(
    i, j, d2, map_matrix,
    index_residue1::Vector{Int},
    index_residue2::Vector{Int},
    discrete::Bool,
)
    (; m) = map_matrix
    ires = index_residue1[i]
    jres = index_residue2[j]
    if discrete
        m[ires, jres] = true
    else
        m[ires, jres] = _contact(m[ires, jres], sqrt(d2))
    end
    return map_matrix
end

function contact_map(
    atoms1::AbstractVector{<:Atom},
    atoms2::AbstractVector{<:Atom};
    dmax::Real=4.0,
    discrete::Bool=true,
    unitcell=nothing,
    positions1::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing,
    positions2::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing,
    parallel::Bool=false,
)
    type = Union{Missing,discrete ? Bool : promote_type(typeof(atoms1[1].x), typeof(atoms2[1].x))}
    residues1 = collect(PDBTools.eachresidue(atoms1))
    index_residue1 = assign_index_residue(atoms1, residues1)
    residues2 = collect(PDBTools.eachresidue(atoms2))
    index_residue2 = assign_index_residue(atoms2, residues2)
    map = zero(ContactMap{type}, length(residues1), length(residues2), dmax, 0, residues1, residues2)
    map_matrix = MapMatrix(map.matrix)
    sys = CellListMap.ParticleSystem(
        xpositions=isnothing(positions1) ? coor.(atoms1) : positions1,
        ypositions=isnothing(positions2) ? coor.(atoms2) : positions2,
        cutoff=dmax,
        unitcell=unitcell,
        output=map_matrix,
        output_name=:map_matrix,
        parallel=parallel,
    )
    CellListMap.map_pairwise!(
        (x, y, i, j, d2, map_matrix) -> update_map_matrix(i, j, d2, map_matrix, index_residue1, index_residue2, discrete),
        sys
    )
    return map
end

@testitem "contact_map" begin
    using PDBTools
    using ShowMethodTesting

    # monomer
    ats = read_pdb(PDBTools.TESTPDB, "protein")
    for parallel in (false, true)
        map = contact_map(ats)
        @test size(map.matrix) == (104, 104)
        @test count(map.matrix) == 1106
        @test parse_show(map) ≈ "ContactMap{Union{Missing, Bool}} of size (104, 104), with threshold 4.0 and gap 0"
        map = contact_map(ats; discrete=false)
        @test sum(skipmissing(map.matrix)) ≈ 2407.2163f0
        @test parse_show(map) ≈ "ContactMap{Union{Missing, Float32}} of size (104, 104), with threshold 4.0 and gap 0"
    end

    # dimer
    ats = read_pdb(PDBTools.DIMERPDB)
    cA = select(ats, "chain A")
    cB = select(ats, "chain B")
    for parallel in (false, true)
        map = contact_map(cA, cB; parallel=parallel)
        @test sum(map.matrix) == 17
        @test map[235, :] == [false, false, true, false, false, false, false, false, false, false, false, false]
        @test count(map[:, 3]) == 3
        @test parse_show(map) ≈ "ContactMap{Union{Missing, Bool}} of size (243, 12), with threshold 4.0 and gap 0"
        map = contact_map(cA, cB; discrete=false, parallel=parallel)
        @test sum(skipmissing(map.matrix)) ≈ 58.00371f0
    end

    # Test with and without periodic boundary conditions
    pdb_pbc = read_pdb(PDBTools.TESTPBC, "protein")
    uc = read_unitcell(PDBTools.TESTPBC)
    pdb_nopbc = read_pdb(PDBTools.TESTNOPBC, "protein")
    for parallel in (false, true)
        pbc = contact_map(pdb_pbc; discrete=true, unitcell=uc, parallel=parallel)
        no_pbc = contact_map(pdb_nopbc; discrete=true, parallel=parallel)
        @test pbc.matrix == no_pbc.matrix
        pbc = contact_map(pdb_pbc; discrete=false, unitcell=uc, parallel=parallel)
        no_pbc = contact_map(pdb_nopbc; discrete=false, parallel=parallel)
        @test all(ismissing(pbc.matrix[i]) | isapprox(pbc.matrix[i], no_pbc.matrix[i]; rtol=1e-2) for i in eachindex(pbc.matrix))
    end

    # Arithmetic operations on contact maps
    pdb = wget("2cpb", "model 1 2")
    m = collect(eachmodel(pdb))
    c1 = contact_map(m[1])
    c2 = contact_map(m[2])
    c3 = c1 + c2
    @test sum(c3.matrix) == 424
    c3 = c2 - c1
    @test sum(skipmissing(c3.matrix)) == -10
    c3 = c1 - c1
    @test sum(skipmissing(c3.matrix)) == 0
    c3 = c1 + c1
    @test sum(c3.matrix) == sum(c1.matrix)
    c2.matrix .= missing
    c3 = c2 - c1
    @test sum(skipmissing(c3.matrix)) == -414
    c3 = c1 - c2
    @test sum(skipmissing(c3.matrix)) == 414

    c1 = contact_map(m[1]; discrete=false)
    c2 = contact_map(m[2]; discrete=false)
    c3 = c1 + c2
    @test sum(skipmissing(c3.matrix)) ≈ 1564.5426f0
    c3 = c2 - c1
    @test sum(skipmissing(c3.matrix)) ≈ -3.5303345f0
    c3 = c1 - c1
    @test sum(skipmissing(c3.matrix)) ≈ 0.0
    c3 = c1 + c1
    @test sum(skipmissing(c3.matrix)) ≈ 2 * sum(skipmissing(c1.matrix))
    c3 = contact_map(m[1]; discrete=false)
    c3.matrix .= missing
    c3 = c3 - c1
    @test all(ismissing, c3.matrix)

    # Errors
    c1 = contact_map(m[1])
    c4 = contact_map(select(m[1], "residue > 1"))
    @test_throws "only defined" c4 - c1
    @test_throws "only defined" c4 + c1
    c1 = contact_map(m[1]; discrete=false)
    c4 = contact_map(select(m[1], "residue > 1"); discrete=false)
    @test_throws "only defined" c4 - c1
    @test_throws "only defined" c4 + c1
    c1 = contact_map(m[1])
    c4 = contact_map(select(m[1]); discrete=false)
    @test_throws "same type" c4 - c1
    @test_throws "same type" c4 + c1

end

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
            @test dpbc ≈ dnopbc rtol = 1e-2
        end
    end
end
