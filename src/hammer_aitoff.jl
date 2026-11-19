export hammer_aitoff_projection, hammer_aitoff_map, HammerAitoffMap

"""
    hammer_aitoff_projection(atoms; center = center_of_mass(atoms)) -> (x, y)

Project the atoms onto a plane using the Hammer-Aitoff equal-area projection.

The coordinates of each atom relative to `center` are converted to a latitude `φ`
and a longitude `λ` (as seen from `center`), which are then mapped to planar
coordinates. The projected points lie inside the ellipse
`(x / 2√2)^2 + (y / √2)^2 ≤ 1`. Atoms located exactly at `center` have no
defined direction and are placed at `(0, 0)`.

The input atoms are not modified.

# Arguments
- `atoms`: a vector of `PDBTools.Atom`.
- `center`: the point from which directions are measured (default: the center of mass).

# Returns
Two vectors with the projected `x` and `y` coordinates of each atom.

# Example

```julia-repl
julia> using PDBTools, LazyArtifacts

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> x, y = hammer_aitoff_projection(atoms);
```
"""
function hammer_aitoff_projection(atoms::AbstractVector{<:Atom}; center = center_of_mass(atoms))
    n = length(atoms)
    x_hammer = zeros(n)
    y_hammer = zeros(n)
    for (i, at) in enumerate(atoms)
        x = at.x - center[1]
        y = at.y - center[2]
        z = at.z - center[3]
        r = sqrt(x^2 + y^2 + z^2)
        # Direction undefined at the center: keep the point at (0, 0)
        iszero(r) && continue
        φ = asin(clamp(z / r, -1, 1)) # latitude,  from -π/2 to π/2
        λ = atan(y, x)                # longitude, from -π to π
        denom = sqrt(1 + cos(φ) * cos(λ / 2))
        x_hammer[i] = 2 * sqrt(2) * cos(φ) * sin(λ / 2) / denom
        y_hammer[i] = sqrt(2) * sin(φ) / denom
    end
    return x_hammer, y_hammer
end

#=
    _binned_statistic_2d(x, y, values; nbins=100, statistic=:mean) -> (xbin, ybin, cbins)

Build a 2D histogram of `values` over the points `(x[k], y[k])`.

The bounding box of the points is divided into `nbins × nbins` equal bins.
Each point is assigned to its bin, and the value of a bin is a statistic
(`:mean`, `:sum` or `:count`) of the `values` of the points that fell inside it.
Bins with no points get `0.0`. Entries of `values` that are `NaN` or `Inf` are ignored.

Returns the bin centers along each axis, `xbin` and `ybin`, and a matrix where
`cbins[j, i]` is the value of the bin with center `(xbin[i], ybin[j])`.
=#
function _binned_statistic_2d(x, y, values; nbins::Integer=100, statistic::Symbol=:mean)
    if !(statistic in (:mean, :sum, :count))
        throw(ArgumentError("statistic must be :mean, :sum or :count, got $statistic"))
    end
    nbins > 0 || throw(ArgumentError("nbins must be positive, got $nbins"))
    if !(length(x) == length(y) == length(values))
        throw(ArgumentError("x, y and values must have the same length, got $(length(x)), $(length(y)) and $(length(values))"))
    end
    isempty(values) && throw(ArgumentError("x, y and values are empty"))

    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    # If all points have the same coordinate, use a unit width to avoid dividing by zero
    xwidth = xmax > xmin ? (xmax - xmin) / nbins : 1.0
    ywidth = ymax > ymin ? (ymax - ymin) / nbins : 1.0

    # Bin centers
    xbin = [xmin + (i - 0.5) * xwidth for i in 1:nbins]
    ybin = [ymin + (j - 0.5) * ywidth for j in 1:nbins]

    # For each bin: sum of the values and number of points inside it
    sums = zeros(nbins, nbins)
    counts = zeros(Int, nbins, nbins)
    for (xk, yk, vk) in zip(x, y, values)
        isfinite(vk) || continue
        # The point with the maximum coordinate falls exactly on the upper
        # edge, so it is placed in the last bin
        i = min(floor(Int, (xk - xmin) / xwidth) + 1, nbins)
        j = min(floor(Int, (yk - ymin) / ywidth) + 1, nbins)
        sums[j, i] += vk
        counts[j, i] += 1
    end

    # Final value of each bin (empty bins stay at 0.0)
    cbins = zeros(nbins, nbins)
    for i in 1:nbins, j in 1:nbins
        counts[j, i] > 0 || continue
        cbins[j, i] = if statistic == :mean
            sums[j, i] / counts[j, i]
        elseif statistic == :sum
            sums[j, i]
        else # :count
            counts[j, i]
        end
    end
    return xbin, ybin, cbins
end

"""
    HammerAitoffMap

Values of a per-atom property on the Hammer-Aitoff projection of a structure, as
returned by [`hammer_aitoff_map`](@ref). Plot it with `heatmap(map)` (requires
`using Plots`).

# Fields
- `x`, `y`: projected coordinates of each atom (see [`hammer_aitoff_projection`](@ref)).
- `chains`: chain of each atom (used to draw the contours of the chains).
- `xbin`, `ybin`: centers of the bins along each axis.
- `values`: matrix where `values[j, i]` is the value of the bin with center `(xbin[i], ybin[j])`.
- `statistic`: how the values inside a bin were combined (`:mean`, `:sum` or `:count`).
"""
struct HammerAitoffMap
    x::Vector{Float64}
    y::Vector{Float64}
    chains::Vector{String}
    xbin::Vector{Float64}
    ybin::Vector{Float64}
    values::Matrix{Float64}
    statistic::Symbol
end

function Base.show(io::IO, ::MIME"text/plain", m::HammerAitoffMap)
    print(io, chomp("""
    HammerAitoffMap of $(length(m.x)) atoms ($(length(unique(m.chains))) chains)
        Bins: $(length(m.xbin)) × $(length(m.ybin)), statistic: :$(m.statistic)
        Range of values: $(extrema(m.values))
    """))
end

"""
    hammer_aitoff_map(atoms, values; nbins=100, statistic=:mean, center=center_of_mass(atoms))

Compute the map of a per-atom property (`values`, one value per atom) on the Hammer-Aitoff
projection of the structure. The atoms are projected with [`hammer_aitoff_projection`](@ref),
the map is divided into `nbins × nbins` bins, and the value of each bin is a statistic of the
values of the atoms that fell inside it. Bins without atoms get `0.0`, and values that are
`NaN` or `Inf` are ignored.

Returns a [`HammerAitoffMap`](@ref), which can be plotted with `heatmap(map)` (requires `using Plots`).

# Keyword arguments
- `nbins`: number of bins along each axis.
- `statistic`: how the values inside a bin are combined:
  - `:mean`: the mean of the values (default), for intensive properties, such as B-factors or distances;
  - `:sum`: the sum of the values, for additive properties, such as atomic contributions to a coordination number;
  - `:count`: the number of atoms in the bin (the values are ignored).
- `center`: the point from which the directions of the atoms are measured.

# Example

```julia-repl
julia> using PDBTools, LazyArtifacts

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> map = hammer_aitoff_map(atoms, mass.(atoms); nbins = 50)
HammerAitoffMap of 1463 atoms (1 chains)
    Bins: 50 × 50, statistic: :mean
    Range of values: (0.0, 32.060001373291016)
```
"""
function hammer_aitoff_map(
    atoms::AbstractVector{<:Atom},
    values::AbstractVector{<:Real};
    nbins::Integer=100,
    statistic::Symbol=:mean,
    center = center_of_mass(atoms),
)
    length(values) == length(atoms) ||
        throw(ArgumentError("values has $(length(values)) elements, but there are $(length(atoms)) atoms"))
    x, y = hammer_aitoff_projection(atoms; center)
    xbin, ybin, cbins = _binned_statistic_2d(x, y, values; nbins, statistic)
    return HammerAitoffMap(x, y, string.(chain.(atoms)), xbin, ybin, cbins, statistic)
end

@testitem "hammer_aitoff_projection" begin
    using PDBTools, LazyArtifacts
    atoms = read_mmcif(PDBTools.test_file("6CO8.cif"), "protein")
    coords_before = [(at.x, at.y, at.z) for at in atoms]
    x, y = hammer_aitoff_projection(atoms)
    @test length(x) == length(y) == length(atoms)
    # the input atoms are not modified
    @test [(at.x, at.y, at.z) for at in atoms] == coords_before
    # every projected point lies inside the Hammer ellipse
    @test all(@. (x / (2 * sqrt(2)))^2 + (y / sqrt(2))^2 <= 1 + 1e-8)
    # projection of known directions (exact coordinates, center at the origin)
    c = [0.0, 0.0, 0.0]
    at = deepcopy(atoms[1])
    at.x, at.y, at.z = 1.0, 0.0, 0.0                # φ = 0, λ = 0
    xp, yp = hammer_aitoff_projection([at]; center = c)
    @test xp[1] ≈ 0.0 atol = 1e-12
    @test yp[1] ≈ 0.0 atol = 1e-12
    at.x, at.y, at.z = 0.0, 0.0, 1.0                # north pole
    xp, yp = hammer_aitoff_projection([at]; center = c)
    @test xp[1] ≈ 0.0 atol = 1e-12
    @test yp[1] ≈ sqrt(2)
end

@testitem "hammer_aitoff_map" begin
    using PDBTools, LazyArtifacts
    using PDBTools: _binned_statistic_2d
    using ShowMethodTesting
    atoms = read_mmcif(PDBTools.test_file("6CO8.cif"), "protein")
    x, y = hammer_aitoff_projection(atoms)
    n = length(x)
    xbin, ybin, c = _binned_statistic_2d(x, y, ones(n); nbins = 20, statistic = :count)
    @test length(xbin) == length(ybin) == 20
    @test size(c) == (20, 20)
    @test sum(c) == n
    # with values = 1, :sum equals :count and every non-empty bin has mean 1
    _, _, s = _binned_statistic_2d(x, y, ones(n); nbins = 20, statistic = :sum)
    @test s == c
    _, _, m = _binned_statistic_2d(x, y, ones(n); nbins = 20, statistic = :mean)
    @test all(m[c .> 0] .== 1.0)
    # non-finite values are ignored
    v = ones(n); v[1] = NaN
    _, _, c2 = _binned_statistic_2d(x, y, v; nbins = 20, statistic = :count)
    @test sum(c2) == n - 1
    # the map contains the projection and the binned values
    map = hammer_aitoff_map(atoms, ones(n); nbins = 20, statistic = :count)
    @test map isa HammerAitoffMap
    @test map.values == c
    @test (map.x, map.y) == (x, y)
    @test map.chains == string.(chain.(atoms))
    @test map.statistic == :count
    @test parse_show(map) ≈ """
        HammerAitoffMap of $n atoms ($(length(unique(chain.(atoms)))) chains)
            Bins: 20 × 20, statistic: :count
            Range of values: $(extrema(c))
        """
    # errors
    @test_throws ArgumentError hammer_aitoff_map(atoms, ones(n); statistic = :median)
    @test_throws ArgumentError hammer_aitoff_map(atoms, ones(n); nbins = 0)
    @test_throws ArgumentError hammer_aitoff_map(atoms, ones(n - 1))
end
