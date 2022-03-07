"""

```
distance(x,y)
```

Computes the minimum distance between two sets of atoms, between an atom and a set of atoms, or simply the distance between two atoms. The input may be a vector of `Atom`s, or the matrices that are output of the `coor` function. Using the matrices may result in greater speed, particularly if the greatest dimension of the matrices is the first one (`xyz_in_cols=true` option, which is the default option).

### Examples

```julia-repl
julia> model = wget("1BSX");

julia> protein = select(model,"protein");

julia> ligand = select(model,"resname T3");

julia> distance(protein,ligand)
2.7775834820937417

julia> distance(protein[1],ligand[3])
36.453551075306784

julia> using BenchmarkTools

julia> @btime distance(\$ligand[1],\$protein)
  10.092 μs (0 allocations: 0 bytes)
4.186781819010877

julia> xprot = coor(protein);

julia> @btime distance(\$ligand[1],\$xprot)
  5.590 μs (0 allocations: 0 bytes)
4.186781819010877

```

"""
@inline distance(x₁, y₁, z₁, x₂, y₂, z₂) = sqrt((x₂ - x₁)^2 + (y₂ - y₁)^2 + (z₂ - z₁)^2)
@inline distance(x, y, z, at::Atom) = distance(x, y, z, at.x, at.y, at.z)
@inline distance(at::Atom, x, y, z) = distance(x, y, z, at)
@inline distance(x::Atom, y::Atom) = distance(x.x, x.y, x.z, y.x, y.y, y.z)

@inline distance_sq(x₁, y₁, z₁, x₂, y₂, z₂) = (x₂ - x₁)^2 + (y₂ - y₁)^2 + (z₂ - z₁)^2
@inline distance_sq(x, y, z, at::Atom) = distance_sq(x, y, z, at.x, at.y, at.z)
@inline distance_sq(at::Atom, x, y, z) = distance_sq(x, y, z, at)
@inline distance_sq(x::Atom, y::Atom) = distance_sq(x.x, x.y, x.z, y.x, y.y, y.z)

"""

```
closest(x,y)
```

Computes the minimum distance between two sets of atoms and returns the index(es) of the atoms and their distance. The input may be a vector of `Atom`s, or the matrices that are output of the `coor` function. Using the matrices results in greater speed, particularly if the greatest dimension of the matrices is the first one (`xyz_in_cols=true` option, which is the default option).

### Examples

```julia-repl
julia> model = wget("1BSX");

julia> protein = select(model,"protein");

julia> ligand = select(model,"resname T3");

julia> closest(ligand,protein)
((43, 3684), 2.7775834820937417)

julia> ligand[43]
    4037   O1      T3     B        2      512  -22.568   81.625    3.159 36.59  1.00     1       -      4041

julia> closest(ligand[1],protein)
(728, 4.186781819010877)

julia> using BenchmarkTools

julia> @btime closest(\$ligand[1],\$protein)
  9.940 μs (0 allocations: 0 bytes)
(728, 4.186781819010877)

julia> xprot = coor(protein);

julia> @btime closest(\$ligand[1],\$xprot)
  5.451 μs (0 allocations: 0 bytes)
(728, 4.186781819010877)


```

"""
function closest(x, y, z, atoms::AbstractVector{Atom})
    d = +Inf
    imin = -1
    for (i, at) in pairs(atoms)
        dᵢ = distance_sq(x, y, z, at)
        if dᵢ < d
            d = dᵢ
            imin = i
        end
    end
    imin, sqrt(d)
end
closest(at::Atom, y::AbstractVector{Atom}) = closest(at.x, at.y, at.z, y)
distance(at::Atom, y::AbstractVector{Atom}) = closest(at, y)[2]

# Two vectors of atoms
function closest(x::AbstractVector{Atom}, y::AbstractVector{Atom})
    lx = length(x)
    ly = length(y)
    if lx >= ly
        outer = x
        inner = y
    else
        outer = y
        inner = x
    end
    d = +Inf
    pair = (-1, -1)
    for (i, atx) in pairs(outer)
        j, dᵢⱼ = closest(atx, inner)
        if dᵢⱼ < d
            d = dᵢⱼ
            pair = (i, j)
        end
    end
    if lx >= ly
        pair, d
    else
        (pair[2], pair[1]), d
    end
end
distance(x::AbstractVector{Atom}, y::AbstractVector{Atom}) = closest(x, y)[2]

# Two residues
closest(r1::Residue, r2::Residue) =
    closest(@view(r1.atoms[r1.range]), @view(r2.atoms[r2.range]))
distance(r1::Residue, r2::Residue) =
    distance(@view(r1.atoms[r1.range]), @view(r2.atoms[r2.range]))

#
# The following versions are thought to be faster because the coordinates
# are provided as matrices
#

# Three coordinates and a matrix of coordinates
function closest(x, y, z, A::Matrix; xyz_in_cols = true)
    n, m = size(A)
    imin = 1
    if xyz_in_cols
        @assert m == 3 "Number of columns of coordinates matrix must be 3 when xyz_in_cols=true"
        d = distance_sq(x, y, z, A[1, 1], A[1, 2], A[1, 3])
        for i = 2:n
            dᵢ = distance_sq(x, y, z, A[i, 1], A[i, 2], A[i, 3])
            if dᵢ < d
                d = dᵢ
                imin = i
            end
        end
    else
        @assert n == 3 "Number of rows of coordinates matrix must be 3 when xyz_in_cols=false"
        d = distance_sq(x, y, z, A[1, 1], A[2, 1], A[3, 1])
        for j = 2:m
            dⱼ = distance_sq(x, y, z, A[1, j], A[2, j], A[3, j])
            if dⱼ < d
                d = dⱼ
                imin = j
            end
        end
    end
    imin, sqrt(d)
end
distance(x, y, z, A::Matrix; xyz_in_cols = true) =
    closest(x, y, z, A, xyz_in_cols = xyz_in_cols)[2]

# One Atom and a matrix of coordinates
closest(at::Atom, A::Matrix; xyz_in_cols = true) =
    closest(at.x, at.y, at.z, A, xyz_in_cols = xyz_in_cols)
distance(at::Atom, A::Matrix; xyz_in_cols = true) =
    distance(at.x, at.y, at.z, A, xyz_in_cols = xyz_in_cols)

#
# Two matrices of coordinates. The smallest matrix is iterated atom by atom,
# and the largest matrix is iterated on the above function, which is loop-vectorized   
#
function closest(A::Matrix, B::Matrix; xyz_in_cols = true)
    nA, mA = size(A)
    nB, mB = size(B)
    pair = (1, 1)
    if xyz_in_cols
        @assert (mA == 3 && mB == 3) "Matrices must have 3 columns if xyz_in_cols=true"
        d = distance(B[1, 1], B[1, 2], B[1, 3], A[1, 1], A[1, 2], A[1, 3])
        if nA > nB
            for i = 2:nB
                j, dᵢⱼ = closest(B[i, 1], B[i, 2], B[i, 3], A, xyz_in_cols = true)
                if dᵢⱼ < d
                    d = dᵢⱼ
                    pair = (i, j)
                end
            end
        else
            for i = 2:nA
                j, dᵢⱼ = closest(A[i, 1], A[i, 2], A[i, 3], B, xyz_in_cols = true)
                if dᵢⱼ < d
                    d = dᵢⱼ
                    pair = (i, j)
                end
            end
        end
    else
        @assert (nA == 3 && nB == 3) "Matrices must have 3 rows if xyz_in_cols=false"
        d = distance(B[1, 1], B[2, 1], B[3, 1], A[1, 1], A[2, 1], A[3, 1])
        if mA > mB
            for j = 2:mB
                i, dᵢⱼ = closest(B[1, j], B[2, j], B[3, j], A, xyz_in_cols = false)
                if dᵢⱼ < d
                    d = dᵢⱼ
                    pair = (i, j)
                end
            end
        else
            for j = 2:mA
                i, dᵢⱼ = closest(A[1, j], A[2, j], A[3, j], B, xyz_in_cols = false)
                if dᵢⱼ < d
                    d = dᵢⱼ
                    pair = (i, j)
                end
            end
        end
    end
    pair, d
end
distance(A::Matrix, B::Matrix; xyz_in_cols = true) =
    closest(A, B, xyz_in_cols = xyz_in_cols)[2]
