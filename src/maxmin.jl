#
# Return the coordinates of the atoms
#
struct MaxMinCoords
    xmin::Vector{Float64}
    xmax::Vector{Float64}
    xlength::Vector{Float64}
end

"""
    maxmin(atoms::Vector{Atom}; selection)

Returns the maximum and minimum coordinates of an atom vector, and the length (maximum minus minimum) in each direction. 

### Example

```julia-repl
julia> protein = wget("1LBD");

julia> maxmin(protein)
 
 Minimum atom coordinates: xmin = [-29.301, 57.178, 45.668]
 Maximum atom coordinates: xmax = [47.147, 99.383, 86.886]
 Length in each direction: xlength = [76.448, 42.205, 41.217999999999996]

```

"""
function maxmin(atoms::AbstractVector{Atom}, selection::String)
    query = parse_query(selection)
    return maxmin(atoms, only=atom -> apply_query(query, atom))
end

function maxmin(atoms::AbstractVector{Atom}; only=all)
    xmin = [+Inf, +Inf, +Inf]
    xmax = [-Inf, -Inf, -Inf]
    for at in atoms
        if only(at)
            xmin .= (min(at.x, xmin[1]), min(at.y, xmin[2]), min(at.z, xmin[3]))
            xmax .= (max(at.x, xmax[1]), max(at.y, xmax[2]), max(at.z, xmax[3]))
        end
    end
    xlength = @. xmax - xmin
    return MaxMinCoords(xmin, xmax, xlength)
end

function Base.show(io::IO, m::MaxMinCoords)
    print(io, chomp("""
    Minimum atom coordinates: xmin = $(m.xmin)
    Maximum atom coordinates: xmax = $(m.xmax)
    Length in each direction: xlength = $(m.xlength)
    """))
end

@testitem "maxmin" begin
    atoms = read_pdb(PDBTools.TESTPDB, "protein")
    m = maxmin(atoms)
    @test m.xmin ≈ [-14.18, -17.561, -15.369]
    @test m.xmax ≈ [18.694, 14.182, 15.909]
    @test m.xlength ≈ [32.873999999999995, 31.743000000000002, 31.278]
end