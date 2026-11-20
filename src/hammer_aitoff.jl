export hammer, hand_made_histogram, hammeraitoff, plot_hammeraitoff

"""
    hammer(atoms; center = center_of_mass(atoms)) -> (x_hammer, y_hammer)

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
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> x, y = hammer(atoms);
```
"""
function hammer(atoms::AbstractVector{<:Atom}; center = center_of_mass(atoms))
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

"""
    hand_made_histogram(x, y, values; nbins=100, statistic=:mean) -> (xbin, ybin, cbins)

Build a 2D histogram of `values` over the points `(x[k], y[k])`.

The bounding box of the points is divided into `nbins × nbins` equal bins.
Each point is assigned to its bin, and the value of a bin is a statistic of
the `values` of the points that fell inside it. Bins with no points get `0.0`.
Entries of `values` that are `NaN` or `Inf` are ignored.

# Arguments
- `x`, `y`: coordinates of the points.
- `values`: one value per point.
- `nbins`: number of bins along each axis.
- `statistic`: how the values inside a bin are combined, one of
  - `:mean`: the mean of the values (default);
  - `:sum`: the sum of the values;
  - `:count`: the number of points in the bin.

# Returns
- `xbin`, `ybin`: the bin centers along each axis.
- `cbins`: a matrix where `cbins[j, i]` is the value of the bin
  with center `(xbin[i], ybin[j])`.
"""
function hand_made_histogram(x, y, values; nbins::Integer=100, statistic::Symbol=:mean)
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
    hammeraitoff(atoms, values; nbins=100, statistic=:mean) -> (xbin, ybin, cbins)

Project `atoms` with [`hammer`](@ref) and build the 2D histogram of `values`
with [`hand_made_histogram`](@ref). The `statistic` can be `:mean`, `:sum`
or `:count`.
"""
function hammeraitoff(atoms::AbstractVector{<:Atom}, values; nbins::Integer=100, statistic::Symbol=:mean)
    length(values) == length(atoms) ||
        throw(ArgumentError("values has $(length(values)) elements, but there are $(length(atoms)) atoms"))
    x_hammer, y_hammer = hammer(atoms)
    return hand_made_histogram(x_hammer, y_hammer, values; nbins, statistic)
end

"""
    plot_hammeraitoff(atoms, values;
                      nbins=100,
                      statistic=:mean,
                      colormap=:tempo,
                      clims=nothing,
                      contours=:none,
                      groupfun=identity,
                      contour_color=:grey30,
                      colorbar_title="",
                      title="",
                      outfile=nothing)

Plot a heat map of `values` (one value per atom) using the Hammer-Aitoff projection.

!!! note
    This function requires `using Plots`.

# Keyword arguments
- `nbins`: number of bins along each axis.
- `statistic`: how the values inside a bin are combined: `:mean`, `:sum`
  or `:count` (see [`hand_made_histogram`](@ref)).
- `colormap`: color map of the heat map.
- `clims`: color limits as a tuple `(min, max)`. If `nothing`, they are
  taken from the data.
- `contours`: contour mode, one of
  - `:none`: no contours;
  - `:chain` (or `:chains`): draw one contour line around each chain;
  - `:group`: draw contours around groups of chains, where the group of a
    chain is given by `groupfun(chainname)`.
- `groupfun`: function mapping a chain name to its group (only used when
  `contours=:group`).
- `contour_color`: color of the contour lines.
- `colorbar_title`, `title`: texts of the color bar and of the plot.
- `outfile`: if given, the figure is saved to this path.

# Returns
The plot object.

# Example

```julia-repl
julia> using PDBTools, Plots

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> plot_hammeraitoff(atoms, beta.(atoms); contours = :chain, colorbar_title = "B-factor")
```
"""
function plot_hammeraitoff end

@testitem "hammer" begin
    using PDBTools, LazyArtifacts
    atoms = read_mmcif(PDBTools.test_file("6CO8.cif"), "protein")
    coords_before = [(at.x, at.y, at.z) for at in atoms]
    x, y = hammer(atoms)
    @test length(x) == length(y) == length(atoms)
    # the input atoms are not modified
    @test [(at.x, at.y, at.z) for at in atoms] == coords_before
    # every projected point lies inside the Hammer ellipse
    @test all(@. (x / (2 * sqrt(2)))^2 + (y / sqrt(2))^2 <= 1 + 1e-8)
    # projection of known directions (exact coordinates, center at the origin)
    c = [0.0, 0.0, 0.0]
    at = deepcopy(atoms[1])
    at.x, at.y, at.z = 1.0, 0.0, 0.0                # φ = 0, λ = 0
    xp, yp = hammer([at]; center = c)
    @test xp[1] ≈ 0.0 atol = 1e-12
    @test yp[1] ≈ 0.0 atol = 1e-12
    at.x, at.y, at.z = 0.0, 0.0, 1.0                # north pole
    xp, yp = hammer([at]; center = c)
    @test xp[1] ≈ 0.0 atol = 1e-12
    @test yp[1] ≈ sqrt(2)
end

@testitem "hand_made_histogram" begin
    using PDBTools, LazyArtifacts
    atoms = read_mmcif(PDBTools.test_file("6CO8.cif"), "protein")
    x, y = hammer(atoms)
    n = length(x)
    xbin, ybin, c = hand_made_histogram(x, y, ones(n); nbins = 20, statistic = :count)
    @test length(xbin) == length(ybin) == 20
    @test size(c) == (20, 20)
    @test sum(c) == n
    # with values = 1, :sum equals :count and every non-empty bin has mean 1
    _, _, s = hand_made_histogram(x, y, ones(n); nbins = 20, statistic = :sum)
    @test s == c
    _, _, m = hand_made_histogram(x, y, ones(n); nbins = 20, statistic = :mean)
    @test all(m[c .> 0] .== 1.0)
    # non-finite values are ignored
    v = ones(n); v[1] = NaN
    _, _, c2 = hand_made_histogram(x, y, v; nbins = 20, statistic = :count)
    @test sum(c2) == n - 1
    # hammeraitoff gives the same result as hammer + hand_made_histogram
    @test hammeraitoff(atoms, ones(n); nbins = 20, statistic = :count)[3] == c
    # errors
    @test_throws ArgumentError hand_made_histogram(x, y, ones(n); statistic = :median)
    @test_throws ArgumentError hand_made_histogram(x, y, ones(n - 1))
    @test_throws ArgumentError hammeraitoff(atoms, ones(n - 1))
end
