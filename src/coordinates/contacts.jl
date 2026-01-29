import CellListMap
using SparseArrays: SparseMatrixCSC, sparse, spzeros

export ContactMap, contact_map, residue_residue_distance

"""
    ContactMap{Bool|Real}

Data structure to store contact maps between residues in a protein structure.
The contact map is a matrix of distances between residues. A contact is defined 
when the distance between any two atoms of the residues is less than a given threshold `dmax`.

If the distance between two residues is greater than `dmax`, there will be no
entry associated with the pair in the sparse array, indicating that there is no 
contact between the residues. The return value of `getindex` will be `missing`.

If the distance is less than `dmax`, the value in the matrix is the inverse of the distance
between the residues, and the return value of `getindex` will be the distance.

The `gap` parameter is used to calculate contacts between residues separated by
a given number of residues. For example, if `gap=3`, the contact map was 
calculated between residues separated by at least 3 residues in the sequence.

# Fields

- `matrix::SparseMatrixCSC{T}`: Sparse matrix of the inverse of the distances between residues. 
- `d::T`: Threshold distance for a contact.
- `gap::Int`: Gap between residues to calculate contacts.

If the contact map was calculated with `discrete=true`, the matrix contains
`Bool` values, where `true` indicates a contact.
On the other hand, if `discrete=false`, the matrix contains the inverse of the distances between
residues. Residue pairs without contacts within the `dmax` will not be stored
in the sparse matrix representation. 

"""
struct ContactMap{T<:Union{Bool,<:Real}}
    matrix::SparseMatrixCSC{T}
    d::Float64
    gap::Int
    residues1::Vector{PDBTools.Residue}
    residues2::Vector{PDBTools.Residue}
end
Base.setindex!(map::ContactMap, value, i, j) = map.matrix[i, j] = value
function Base.getindex(map::ContactMap{T}, i, j) where {T<:Real}
    value = map.matrix[i, j]
    if value == zero(T) 
        missing
    else
        inv(value)
    end
end
Base.getindex(map::ContactMap{Bool}, i, j) = map.matrix[i,j]

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

function +(c1::ContactMap{Bool}, c2::ContactMap{Bool})
    _err_sum_cmap(c1, c2)
    c3_matrix = c1.matrix .| c2.matrix
    return ContactMap(c3_matrix, c1.d, c1.gap, c1.residues1, c1.residues2)
end
function -(c1::ContactMap{Bool}, c2::ContactMap{Bool})
    _err_sum_cmap(c1, c2)
    c3_matrix = spzeros(Int, size(c1.matrix))
    for i in eachindex(c1.matrix, c2.matrix, c3_matrix)
        if (c1.matrix[i] === true) & (!(c2.matrix[i] === true))
            c3_matrix[i] = 1
        elseif (!(c1.matrix[i] === true)) & (c2.matrix[i] === true)
            c3_matrix[i] = -1
        end
    end
    return ContactMap(c3_matrix, c1.d, c1.gap, c1.residues1, c1.residues2)
end

function +(c1::ContactMap{T1<:AbstractFloat}, c2::ContactMap{T2<:AbstractFloat}) where {T1,T2}
    _err_sum_cmap(c1, c2)
    T = promote_type(T1,T2)
    c1_i, c1_j, c1_v = findnz(c1)
    c2_i, c2_j, c2_v = findnz(c2)
    c_sum = sparse(
        vcat(c1_i, c2_i),
        vcat(c1_j, c2_j),
        vcat(c1_v, c2_v),
        length(c1.residues1), length(c1.residues2), 
        (x,y) -> _combine(T,x,y)
    )
    return ContactMap(c_sum, c1.d, c1.gap, c1.residues1, c1.residues2)
end
voltar

function -(c1::ContactMap{<:AbstractFloat}, c2::ContactMap{<:AbstractFloat})
    _err_sum_cmap(c1, c2)
    c_sub = @. inv(inv(c1) - inv(c2)) 
    return ContactMap(c_sub, c1.d, c1.gap, c1.residues1, c1.residues2)
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

julia> map = contact_map(cA) # contact map of chain A
ContactMap{Bool} of size (243, 243), with threshold 4.0 and gap 0 

julia> map[95,95] # fetch the existence of a contact
true

julia> # using Plots; heatmap(map); # uncomment to plot the contact map
```

## Contact map between residues in two different structures

```jldoctest
julia> using PDBTools

julia> ats = read_pdb(PDBTools.DIMERPDB);

julia> cA = select(ats, "chain A");

julia> cB = select(ats, "chain B");

julia> map = contact_map(cA, cB) # contact map between chains A and B
ContactMap{Bool} of size (243, 12), with threshold 4.0 and gap 0 

julia> map[80, 5] # fetch existence of contact 
true

julia> # using Plots; heatmap(map); # uncomment plot the contact map
```

"""
function contact_map end

#
# Fast contact map calculation for big structures using cell lists.
#
struct MapMatrix{T}
    col_ind::Vector{Int}
    row_ind::Vector{Int}
    value::Vector{T}
end

function _contact!(m::MapMatrix, ires, jres; symmetric=false)
    push!(m.col_ind, ires)
    push!(m.row_ind, jres)
    push!(m.value, true)
    if symmetric
        push!(m.row_ind, ires)
        push!(m.col_ind, jres)
        push!(m.value, true)
    end
    return m
end

function _contact!(m::MapMatrix, ires, jres, d; symmetric=false)
    push!(m.col_ind, ires)
    push!(m.row_ind, jres)
    push!(m.value, d)
    if symmetric
        push!(m.row_ind, ires)
        push!(m.col_ind, jres)
        push!(m.value, d)
    end
    return m
end

CellListMap.copy_output(x::MapMatrix) = MapMatrix(copy(x.col_ind), copy(x.row_ind), copy(x.value))
function CellListMap.reducer(x::MapMatrix, y::MapMatrix)
    append!(x.row_ind, y.row_ind)
    append!(x.col_ind, y.col_ind)
    append!(x.value, y.value)
    return x
end
function CellListMap.reset_output!(x::MapMatrix)   
    empty!(x.row_ind)
    empty!(x.col_ind)
    empty!(x.value)
    return x
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

function update_map_matrix(
    i, j, d2, map_matrix,
    index_residue::Vector{Int},
    gap::Real,
    discrete::Bool
)
    ires = index_residue[i]
    jres = index_residue[j]
    d_gap = abs(ires - jres)
    if d_gap >= gap
        if discrete
            _contact!(map_matrix, ires, jres; symmetric=true)
        else
            if d_gap == 0 # same-residue distance is zero
                _contact!(map_matrix, ires, jres, inv(zero(d2)); symmetric=true)
            else
                _contact!(map_matrix, ires, jres, inv(sqrt(d2)); symmetric=true)
            end
        end
    end
    return map_matrix
end

# 
# Function to combine elements of same indices in the sparse matrix
# final construction: the stored value is the *inverse* of the distance,
# such that 0-values are equivalent to non-stored values in the sparse
# matrix representation.
#
_combine(::Type{Bool}, x, y) = x | y
_combine(::Type{<:Real}, x, y) = max(x,y)

function contact_map(
    atoms1::AbstractVector{<:PDBTools.Atom};
    dmax::Real=4.0,
    gap::Int=0,
    discrete::Bool=true,
    unitcell=nothing,
    positions::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing,
    parallel::Bool=false,
)
    T = discrete ? Bool : typeof(atoms1[1].x)
    residues = collect(PDBTools.eachresidue(atoms1))
    index_residue = assign_index_residue(atoms1, residues)
    map_matrix = MapMatrix(Int[], Int[], T[])
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
    map = ContactMap{T}(
        sparse(
            map_matrix.col_ind, map_matrix.row_ind, map_matrix.value, 
            length(residues), length(residues), 
            (x,y) -> _combine(T,x,y)
        ), 
        dmax, 
        gap, 
        residues, 
        residues
    )
    return map
end

function update_map_matrix(
    i, j, d2, map_matrix,
    index_residue1::Vector{Int},
    index_residue2::Vector{Int},
    discrete::Bool,
)
    ires = index_residue1[i]
    jres = index_residue2[j]
    if discrete
        _contact!(map_matrix, ires, jres)
    else
        _contact!(map_matrix, ires, jres, inv(sqrt(d2)))
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
    T = discrete ? Bool : promote_type(typeof(atoms1[1].x), typeof(atoms2[1].x))
    residues1 = collect(PDBTools.eachresidue(atoms1))
    index_residue1 = assign_index_residue(atoms1, residues1)
    residues2 = collect(PDBTools.eachresidue(atoms2))
    index_residue2 = assign_index_residue(atoms2, residues2)
    map_matrix = MapMatrix(Int[], Int[], T[])
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
    map = ContactMap{T}(
        sparse(
            map_matrix.col_ind, map_matrix.row_ind, map_matrix.value, 
            length(residues1), length(residues2), 
            (x,y) -> _combine(T,x,y)
        ), 
        dmax, 
        0, 
        residues1, 
        residues2,
    )
    return map
end

@testitem "contact_map" begin
    using PDBTools
    using ShowMethodTesting
    using SparseArrays

    # monomer
    ats = read_pdb(PDBTools.TESTPDB, "protein")
    for parallel in (false, true)
        map = contact_map(ats)
        @test size(map.matrix) == (104, 104)
        @test count(map.matrix) == 1106
        @test parse_show(map) ≈ "ContactMap{Bool} of size (104, 104), with threshold 4.0 and gap 0"
        map = contact_map(ats; discrete=false)
        @test sum(inv.(nonzeros(map.matrix))) ≈ 2407.2163f0
        @test parse_show(map) ≈ "ContactMap{Float32} of size (104, 104), with threshold 4.0 and gap 0"
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
        @test parse_show(map) ≈ "ContactMap{Bool} of size (243, 12), with threshold 4.0 and gap 0"
        map = contact_map(cA, cB; discrete=false, parallel=parallel)
        @test sum(inv.(nonzeros(map.matrix))) ≈ 58.00371f0
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
        @test all(isapprox(pbc.matrix[i], no_pbc.matrix[i]; rtol=1e-2) for i in eachindex(pbc.matrix))
    end

    # Arithmetic operations on contact maps
    pdb = wget("2cpb", "model 1 2")
    m = collect(eachmodel(pdb))
    c1 = contact_map(m[1])
    c2 = contact_map(m[2])
    c3 = c1 + c2
    @test sum(c3.matrix) == 424
    c3 = c2 - c1
    @test sum(c3.matrix) == -10
    c3 = c1 - c1
    @test sum(c3.matrix) == 0
    c3 = c1 + c1
    @test sum(c3.matrix) == sum(c1.matrix)
    droptol!(c2.matrix, 10)
    c3 = c2 - c1
    @test sum(c3.matrix) == -414
    c3 = c1 - c2
    @test sum(c3.matrix) == 414

    c1 = contact_map(m[1]; discrete=false)
    c2 = contact_map(m[2]; discrete=false)
    voltar
    c3 = c1 + c2
    @test sum(c3.matrix) ≈ 1564.5426f0
    c3 = c2 - c1
    @test sum(c3.matrix) ≈ -3.5303345f0
    c3 = c1 - c1
    @test sum(c3.matrix) ≈ 0.0
    c3 = c1 + c1
    @test sum(c3.matrix) ≈ 2 * sum(c1.matrix)
    c3 = contact_map(m[1]; discrete=false)
    droptol!(c3.matrix, 10)
    c3 = c3 - c1
    @test nnz(c3.matrix) == 0  

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
