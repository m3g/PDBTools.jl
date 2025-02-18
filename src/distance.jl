"""
    distance(x,y)

Computes the minimum distance between two sets of atoms, between an atom and a set of atoms, or simply 
the distance between two atoms, or from the coordinates or sets of coordinates. 

The input may be an `Atom` vector of `Atom`s, or a 3D vector, or a vector of 3D vector of coordinates, for examples
as output by the `coor` function.

### Examples

```julia-repl
julia> model = wget("1BSX");

julia> protein = select(model,"protein");

julia> ligand = select(model,"resname T3");

julia> distance(protein,ligand)
2.7775834820937417

julia> distance(protein[1],ligand[3])
36.453551075306784

julia> distance(coor(ligand),protein)
2.7775834820937417

```

"""
distance(x::Atom, y::Atom) = norm(coor(x) - coor(y))
distance(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) = norm(x - y)
distance(x::Atom, y::AbstractVector{<:Real}) = norm(coor(x) - y)
distance(x::AbstractVector{<:Real}, y::Atom) = norm(x - coor(y))
distance(x, y) = last(closest(x,y))

"""
    closest(x,y)

Computes the minimum distance between two sets of atoms and returns the indices of the atoms 
and their distance. Both vector of atoms or vectors of coordinates can be used as input.

### Examples

```julia-repl
julia> model = wget("1BSX");

julia> protein = select(model,"protein");

julia> ligand = select(model,"resname T3");

julia> closest(ligand,protein)
(43, 3684, 2.7775834820937417)

julia> ligand[43]
    4037   O1      T3     B        2      512  -22.568   81.625    3.159 36.59  1.00     1       -      4041

julia> closest(ligand[43],protein)
(1, 3684, 2.7775834820937417)

julia> x = coor(protein)
3994-element Vector{SVector{3, Float64}}:
 [52.884, 24.022, 35.587]
 [52.916, 24.598, 36.993]
 ⋮
 [-46.887, 86.925, 13.235]
 [-47.164, 83.593, 15.25]

julia> closest(ligand,x)
(43, 3684, 2.7775834820937417)

```

"""
function closest end

#
# Wrap individual atoms or coordinates in a SVector to dispatch to the _closest function
#
# Distances involving single atoms and single coordinates 
closest(x::Atom, y::Atom) = _closest(SVector{1}(x), SVector{1}(y))
closest(x::Atom, y::AbstractVector{<:Real}) = _closest(SVector{1}(x), SVector{1}((y,)))
closest(x::AbstractVector{<:Real}, y::Atom) = _closest(SVector{1}((x,)), SVector{1}(y))
closest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) = _closest(SVector{1}((x,)), SVector{1}((y,)))
# Distances involving single atoms and vectors of atoms or coordinates
closest(x::Atom, y::AbstractVector{<:Atom}) = _closest(SVector{1}(x), y)
closest(x::Atom, y::AbstractVector{<:AbstractVector{<:Real}}) = _closest(SVector{1}(x), y)
closest(x::AbstractVector{<:Atom}, y::Atom) = _closest(x, SVector{1}(y))
closest(x::AbstractVector{<:AbstractVector{<:Real}},y::Atom) = _closest(x, SVector{1}(y))
# Distances involving single coordinates and vectors of atoms or coordinates
closest(x::AbstractVector{<:Real}, y::AbstractVector{<:Atom}) = _closest(SVector{1}((x,)), y)
closest(x::AbstractVector{<:Atom}, y::AbstractVector{<:Real}) = _closest(x, SVector{1}((y,)))
closest(x::AbstractVector{<:Real}, y::AbstractVector{<:AbstractVector{<:Real}}) = _closest(SVector{1}((x,)), y)
closest(x::AbstractVector{<:AbstractVector{<:Real}}, y::AbstractVector{<:Real}) = _closest(x, SVector{1}((y,)))
# Distance between vectors of atoms / vectors of coordinates
closest(x::AbstractVector{<:Atom}, y::AbstractVector{<:Atom}) = _closest(x, y)
closest(x::AbstractVector{<:AbstractVector{<:Real}}, y::AbstractVector{<:AbstractVector{<:Real}}) = _closest(x, y)
# Distances involving vectors of atoms and vector of coordinates
closest(x::AbstractVector{<:Atom}, y::AbstractVector{<:AbstractVector{<:Real}}) = _closest(x, y)
closest(x::AbstractVector{<:AbstractVector{<:Real}}, y::AbstractVector{<:Atom}) = _closest(x, y)
# Distances involving for residues
closest(x::Residue, y::Residue) = _closest(x.atoms[x.range], y.atoms[y.range])
closest(x::Residue, y::Atom) = _closest(x.atoms[x.range], SVector{1}(y))
closest(x::Atom, y::Residue) = _closest(SVector{1}(x), y.atoms[y.range])
closest(x::Residue, y::AbstractVector{<:Real}) = _closest(x.atoms[x.range], SVector{1}((y,)))
closest(x::AbstractVector{<:Real}, y::Residue) = _closest(SVector{1}((x,)), y.atoms[y.range])
closest(x::Residue, y::AbstractVector{<:SVector}) = _closest(x.atoms[x.range], y)
closest(x::AbstractVector{<:SVector}, y::Residue) = _closest(x, y.atoms[y.range])

_float_type(::AbstractVector{<:AbstractVector{T1}}) where {T1<:Real} = T1
_float_type(atoms::AbstractVector{<:Atom}) = eltype(coor(first(atoms)))
_coor(x::AbstractVector{<:Real}) = x
_coor(x::Atom) = coor(x)

function _closest(
    x::AbstractVector{T1}, 
    y::AbstractVector{T2}, 
) where {T1,T2} 
    imin = -1
    jmin = -1
    dmin = typemax(promote_type(_float_type(x), _float_type(y)))
    for (i, xatom) in pairs(x)
        for (j, yatom) in pairs(y)
            d = norm(_coor(xatom) - _coor(yatom))
            if d < dmin
                imin = i
                jmin = j
                dmin = d
            end
        end
    end
    return imin, jmin, dmin
end

@testitem "distance/closest" begin
    using PDBTools
    atoms = read_pdb(PDBTools.TESTPDB)
    s = select(atoms, "residue = 3")
    s2 = select(atoms, "residue = 5")
    @test distance(s, s2) ≈ 3.6750402718881863 atol=1e-3
    x1 = coor(s)
    x2 = coor(s2)
    @test distance(x1, x2) ≈ 3.6750402718881863 atol=1e-3
    residues = collect(eachresidue(atoms))
    @test distance(residues[3], residues[5]) ≈ 3.6750402718881863 atol=1e-3

    #
    # Dispatch of closest and distance functions 
    #
    r1 = select(atoms, "residue = 3")
    r2 = select(atoms, "residue = 5")

    @test all(isapprox.(closest(r1, r2),(11, 2, 3.6750402718881863); atol=1e-3))
    @test all(isapprox.(closest(r1, coor(r2)),(11, 2, 3.6750402718881863); atol=1e-3))
    @test all(isapprox.(closest(coor(r1), coor(r2)),(11, 2, 3.6750402718881863); atol=1e-3))
    @test all(isapprox.(closest(r1, coor(r2[1])), (11, 1, 3.9481035953986816); atol=1e-3))
    @test all(isapprox.(closest(coor(r1[1]), r2[1]), (1, 1, 5.121629623469468); atol=1e-3))
    @test all(isapprox.(closest(coor(r1), r2[1]),(11, 1, 3.9481035953986816); atol=1e-3))
    @test all(isapprox.(closest(coor(r1), coor(r2[1])),(11, 1, 3.9481035953986816); atol=1e-3))
    @test all(isapprox.(closest(r1[1], coor(r2)),(1, 2, 5.121218702613667); atol=1e-3))
    @test all(isapprox.(closest(coor(r1[1]), coor(r2)),(1, 2, 5.121218702613667); atol=1e-3))
    @test all(isapprox.(closest(coor(r1[1]), r2),(1, 2, 5.121218702613667); atol=1e-3))
    @test all(isapprox.(closest(coor(r1[1]), coor(r2[2])),(1, 1, 5.121218702613667); atol=1e-3))
    @test all(isapprox.(closest(atoms[1], atoms[2]), (1, 1, 0.9994303377424563); atol=1e-3))
    @test all(isapprox.(closest(atoms[1], coor(atoms[2])),(1, 1, 0.9994303377424563); atol=1e-3))
    @test all(isapprox.(closest(coor(atoms[1]), atoms[2]),(1, 1, 0.9994303377424563); atol=1e-3))
    @test all(isapprox.(closest(coor(atoms[1]), coor(atoms[2])), (1, 1, 0.9994303377424563); atol=1e-3))

    @test all(isapprox.(closest(r1[1], coor(r2[2])),(1, 1, 5.121218702613667); atol=1e-3))
    @test all(isapprox.(closest(coor(r1[1]), r2[2]),(1, 1, 5.121218702613667); atol=1e-3))
    @test all(isapprox.(closest(r1[1], r2[2]),(1, 1, 5.121218702613667); atol=1e-3))
    @test all(isapprox.(closest(r1[1], r2), (1, 2, 5.121218702613667); atol=1e-3))
    @test all(isapprox.(closest(r2, r1[1]), (2, 1, 5.121218702613667); atol=1e-3))

    @test distance(r1, r2) ≈ 3.6750402718881863 atol=1e-3
    @test distance(coor(r1), r2) ≈ 3.6750402718881863 atol=1e-3
    @test distance(r1, coor(r2)) ≈ 3.6750402718881863 atol=1e-3
    @test distance(coor(r1), coor(r2)) ≈ 3.6750402718881863 atol=1e-3

    @test distance(r1[1], coor(r2)) ≈ 5.121218702613667 atol=1e-3
    @test distance(coor(r1[1]), coor(r2)) ≈ 5.121218702613667 atol=1e-3
    @test distance(coor(r1[1]), r2) ≈ 5.121218702613667 atol=1e-3
    @test distance(coor(r1[1]), coor(r2[2])) ≈ 5.121218702613667 atol=1e-3

    @test distance(r1[1], coor(r2[2])) ≈ 5.121218702613667 atol=1e-3
    @test distance(coor(r1[1]), r2[2]) ≈ 5.121218702613667 atol=1e-3
    @test distance(r1[1], r2[2]) ≈ 5.121218702613667 atol=1e-3

    r = collect(eachresidue(atoms))
    @test all(isapprox.(closest(r[1], [0.0, 0.0, 0.0]), (12, 1, 16.545482827648158); atol=1e-3))
    @test all(isapprox.(closest([0.0, 0.0, 0.0], r[1]), (1, 12, 16.545482827648158); atol=1e-3))
    @test all(isapprox.(closest(Atom(), r[1]), (1, 12, 16.545482827648158); atol=1e-3))
    @test all(isapprox.(closest(r[1], Atom()), (12, 1, 16.545482827648158); atol=1e-3))
    @test all(isapprox.(closest([coor(Atom())], r[1]), (1, 12, 16.545482827648158); atol=1e-3))
    @test all(isapprox.(closest(r[1], [coor(Atom())]), (12, 1, 16.545482827648158); atol=1e-3))

end

