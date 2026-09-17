#=

Cavity exclusion for sasa_particles
====================================

Optional, opt-in refinement of the Shrake-Rupley-style dot algorithm implemented in
sasa.jl. That algorithm tests, for each atom, whether each of its surface dots is
covered by any *single* neighboring atom's inflated sphere -- a purely local, pairwise
test with no notion of whether a dot that survives that test is actually connected to
bulk solvent, or sits inside a sealed interior cavity. This is the same category of
algorithm used by, e.g., GROMACS's `gmx sasa` (an Eisenhaber et al. 1995 "double cubic
lattice" variant of Shrake-Rupley) and VMD's `measure sasa`, and for typical protein
structures the resulting error is small (a fraction of a percent to a couple of percent
of the total SASA), because most proteins don't have much topologically sealed interior
void space at all.

Some reference programs -- notably SurfaceRacer (Tsodikov, Record & Sergeev, 2002),
used to compute the ASA underlying several of the transfer-model alpha values in
Record.jl -- instead compute "outside" (solvent-connected) ASA as topologically
distinct from interior-cavity ASA, and explicitly exclude the latter. Reproducing that
specific convention is the sole purpose of this file: it is a compatibility mode, not a
general correctness fix (see the discussion in the MTRecord/richards.jl tests for the
empirical evidence on when the two conventions agree and when they don't).

Implementation: connectivity is determined directly on the dots the pairwise test
already found exposed -- not on a separate discretized voxel grid. An earlier version
of this file used a voxel-grid flood fill (inflate each atom onto a 3D grid, flood fill
from a padded, guaranteed-atom-free boundary, and test each dot's grid cell for
reachability); it is not used anymore because it was both less accurate and far more
sensitive to its own tuning parameters than the approach below. On the 3CNA
tetramer/dimer benchmark used throughout this file's tests, the voxel method recovered
only about a fifth of the gap between the uncorrected default and SurfaceRacer's own
number, and its result changed non-monotonically (and sometimes drastically) with grid
spacing and search radius. That is an inherent property of trying to resolve genuinely
narrow-but-real solvent channels and genuinely sealed cavities with a single fixed grid
resolution: whatever window is wide enough to bridge a real channel that discretization
happened to wall off is, for a while, also wide enough to reach across into a real
cavity nearby, and wide enough again eventually washes out cavity detection entirely.

The approach here has no such tension because it never discretizes space at all:

1. Collect the (already atom-atom-pairwise-tested) exposed dots of every atom into one
   point cloud, in absolute coordinates.
2. Connect two dots (via union-find) whenever they are within `cavity_dot_cutoff` of
   each other -- i.e., reconstruct the adjacency of the actual computed molecular
   surface from the dot cloud itself, rather than from a resampled grid.
3. Seed the "exterior" component with the dots most extreme along each Cartesian axis
   (the six that individually maximize/minimize x, y, or z): each of these is, by
   construction, on the true convex, solvent-exposed exterior, and using six independent
   seeds (rather than one) guards against any single one of them landing in a small
   disconnected sliver by coincidence.
4. Any exposed dot whose connected component does not contain one of these seeds is
   judged to be in a sealed interior cavity and is excluded.

Because step 2 operates on the real dot positions (whose local spacing already reflects
each atom's own radius and `n_dots`, via `generate_dots`), the cutoff only needs to be a
small multiple of that spacing to bridge adjacent dots -- with a wide plateau of
essentially indistinguishable results (0.6-1.3 Å tested on the 3CNA benchmark, for the
package's default `n_dots=512`), rather than the voxel method's narrow, non-monotonic
sweet spot.

This module only ever *removes* dots that the pairwise test already marked exposed; it
never adds any. It has no notion of periodic boundary conditions.

=#

#=
    _cavity_seed_indices(positions)

Returns up to 6 indices into `positions` (a `Vector{SVector{3,Float32}}`): the dots
that individually maximize or minimize each of x, y, z. Each is guaranteed to lie on
the true convex exterior of the structure, and is used as an independent connectivity
seed for the "definitely exterior" component (see this file's header comment for why
more than one seed is used).
=#
function _cavity_seed_indices(positions::Vector{SVector{3,Float32}})
    n = length(positions)
    n == 0 && return Int[]
    seeds = Set{Int}()
    for dim in 1:3
        push!(seeds, argmax(i -> positions[i][dim], 1:n))
        push!(seeds, argmin(i -> positions[i][dim], 1:n))
    end
    return collect(seeds)
end

#=
    _cavity_dot_spacing(atom_type, atom_radius_from_type, probe_radius, atoms, n_dots)

Estimates a representative nearest-neighbor spacing between dots on a
`generate_dots`-covered sphere, used to set a sensible default `cavity_dot_cutoff`: the
average area per dot on a sphere of the largest inflated atom radius present is
`4*pi*radius^2 / n_dots`, and the spacing between dots at that density is the square
root of that area.
=#
function _cavity_dot_spacing(atoms, atom_type::Function, atom_radius_from_type::Function, probe_radius::Real, n_dots::Integer)
    rmax = maximum(Float32(atom_radius_from_type(atom_type(at))) + Float32(probe_radius) for at in atoms)
    return sqrt(4.0f0 * Float32(pi) * rmax^2 / n_dots)
end

#=
    _collect_exposed_dot_positions(surface_dots, atoms, dot_cache, atom_type)

Collects the (already atom-atom-pairwise-tested) exposed dots of every atom into one
point cloud, in absolute coordinates, together with back-references to the owning
(atom, dot) index. Split out of `exclude_cavity_dots!` as its own function -- rather
than inlined as its first phase -- purely as a compiler-performance function barrier:
kept together with the (much larger, closure-heavy) connectivity phase below in a
single function, this collection loop measured ~15x slower and allocated ~100x more
than it does on its own, apparently because the combined function was too large for
Julia/LLVM's escape analysis to prove any of the closures in the connectivity phase
non-escaping. Splitting the two phases at a function boundary restores that.
=#
function _collect_exposed_dot_positions(surface_dots, atoms, dot_cache, atom_type::Function)
    positions = SVector{3,Float32}[]
    owner_atom = Int[]
    owner_dot = Int[]
    for i in eachindex(atoms)
        at = atoms[i]
        atom_pos = SVector{3,Float32}(at.x, at.y, at.z)
        dc = dot_cache[atom_type(at)]
        exposed_i = surface_dots[i].exposed
        for idot in eachindex(exposed_i)
            exposed_i[idot] || continue
            push!(positions, atom_pos + SVector{3,Float32}(dc.x[idot], dc.y[idot], dc.z[idot]))
            push!(owner_atom, i)
            push!(owner_dot, idot)
        end
    end
    return positions, owner_atom, owner_dot
end

#=
    _CavityPairs

Per-thread output accumulator for the `pairwise!` neighbor search in
`_exclude_disconnected_dots!`: just the pairs of exposed-dot indices found within
`cutoff` of each other. Implements the `CellListMap` output protocol
(`copy_output`/`reset_output!`/`reducer`, the same one `AtomDotMatrix` in sasa.jl uses)
so the search itself can run multi-threaded; the union-find merge of the collected
pairs is then done serially in `_exclude_disconnected_dots!`, since union-find over a
shared `parent` array is not safely parallelizable across threads the way independent
per-pair accumulation is.
=#
struct _CavityPairs
    pairs::Vector{Tuple{Int,Int}}
end
CellListMap.copy_output(p::_CavityPairs) = _CavityPairs(copy(p.pairs))
function CellListMap.reset_output!(p::_CavityPairs)
    empty!(p.pairs)
    return p
end
function CellListMap.reducer(x::_CavityPairs, y::_CavityPairs)
    append!(x.pairs, y.pairs)
    return x
end
_collect_cavity_pair!(pair, pairs::_CavityPairs) = (push!(pairs.pairs, (pair.i, pair.j)); pairs)

#=
    _find_cavity_pairs(positions, cutoff; parallel=true)

Returns every pair of indices into `positions` that are within `cutoff` of each other,
via `CellListMap.pairwise!` -- the same infrastructure `sasa_particles` uses for the
atom-atom dot-occlusion test -- rather than a hand-rolled `Dict`-based spatial hash: on
the 3CNA benchmark this cut the search from ~11.5 ms to ~1.7-4 ms (parallel/serial) for
~23k exposed dots, since CellListMap's cell list avoids per-lookup tuple hashing and
(when `parallel=true`) threads the search -- worthwhile here since this runs once per
frame of an MD trajectory. Kept as its own function (rather than inlined into
`_exclude_disconnected_dots!` alongside the union-find below) for the same
compiler-performance reason documented on `_collect_exposed_dot_positions`: combined
into one function, the two allocated ~65x more and ran ~1.4x slower, apparently for the
same reason -- too large a function for escape analysis to prove the union-find
closures non-escaping.
=#
function _find_cavity_pairs(positions::Vector{SVector{3,Float32}}, cutoff::Float32; parallel::Bool=true)
    system = ParticleSystem(
        xpositions=positions,
        unitcell=nothing,
        cutoff=cutoff,
        output=_CavityPairs(Tuple{Int,Int}[]),
        output_name=:pairs,
        parallel=parallel,
    )
    pairwise!((pair, pairs) -> _collect_cavity_pair!(pair, pairs), system)
    return system.pairs.pairs
end

#=
    _connected_components(n, pairs)

Plain union-find over `1:n`, unioning each `(i, j)` in `pairs`: returns a `Vector{Int}`
mapping each index to its (fully path-compressed) component root. Kept as its own
function for the same compiler-performance reason documented on `_find_cavity_pairs`.
=#
function _connected_components(n::Integer, pairs::Vector{Tuple{Int,Int}})
    parent = collect(1:n)
    function _find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end
    function _unite!(x, y)
        rx, ry = _find(x), _find(y)
        rx != ry && (parent[rx] = ry)
        return nothing
    end
    for (i, j) in pairs
        _unite!(i, j)
    end
    for i in 1:n
        parent[i] = _find(i)
    end
    return parent
end

#=
    _exclude_disconnected_dots!(surface_dots, positions, owner_atom, owner_dot, cutoff; parallel=true)

Given the exposed-dot point cloud collected by `_collect_exposed_dot_positions`, clears
every dot (via `surface_dots[owner_atom[i]].exposed[owner_dot[i]] = false`) that is not
connected -- via a chain of other exposed dots no farther than `cutoff` apart -- to one
of the structure's extreme (guaranteed-exterior) dots. Mutates `surface_dots` in place
and also returns it. Kept as its own function for the same compiler-performance reason
documented on `_collect_exposed_dot_positions`.
=#
function _exclude_disconnected_dots!(surface_dots, positions, owner_atom, owner_dot, cutoff::Float32; parallel::Bool=true)
    n = length(positions)
    n == 0 && return surface_dots

    pairs = _find_cavity_pairs(positions, cutoff; parallel)
    component = _connected_components(n, pairs)

    exterior_roots = Set(component[i] for i in _cavity_seed_indices(positions))

    for i in 1:n
        if component[i] ∉ exterior_roots
            surface_dots[owner_atom[i]].exposed[owner_dot[i]] = false
        end
    end
    return surface_dots
end

#=
    exclude_cavity_dots!(surface_dots, atoms, dot_cache, atom_type, atom_radius_from_type, probe_radius; n_dots, cavity_dot_cutoff=nothing)

Given the `surface_dots` output of the pairwise dot-occlusion test already run by
`_compute_sasa_particles` (an object indexable as `surface_dots[i].exposed`, an
`AbstractVector{Bool}` of per-dot exposure for atom `i`), clears every exposed dot that
is not connected -- via a chain of other exposed dots no farther than
`cavity_dot_cutoff` apart -- to one of the structure's extreme (guaranteed-exterior)
dots. Mutates `surface_dots` in place and also returns it.

`cavity_dot_cutoff` defaults to twice the estimated nearest-neighbor dot spacing (see
`_cavity_dot_spacing`); results are insensitive to the exact multiple over a wide range
(see this file's header comment), but a value much smaller than the true dot spacing
will spuriously disconnect ordinary exposed surface, and a value much larger will start
bridging across genuine cavities.
=#
function exclude_cavity_dots!(
    surface_dots,
    atoms,
    dot_cache,
    atom_type::Function,
    atom_radius_from_type::Function,
    probe_radius::Real;
    n_dots::Integer,
    cavity_dot_cutoff::Union{Nothing,Real}=nothing,
    parallel::Bool=true,
)
    cutoff = Float32(isnothing(cavity_dot_cutoff) ?
        2 * _cavity_dot_spacing(atoms, atom_type, atom_radius_from_type, probe_radius, n_dots) :
        cavity_dot_cutoff
    )
    cutoff > 0 || throw(ArgumentError("cavity_dot_cutoff must be positive, got $cutoff."))

    positions, owner_atom, owner_dot = _collect_exposed_dot_positions(surface_dots, atoms, dot_cache, atom_type)
    return _exclude_disconnected_dots!(surface_dots, positions, owner_atom, owner_dot, cutoff; parallel)
end

@testitem "exclude_cavities: opt-in, off by default" begin
    using PDBTools
    prot = select(read_pdb(PDBTools.TESTPDB), "protein")
    s_default = sasa_particles(prot)
    s_explicit_off = sasa_particles(prot; exclude_cavities=false)
    @test sasa(s_default) == sasa(s_explicit_off)

    # A single, small, ~104-residue monomeric domain has essentially no topologically
    # sealed interior void space (nothing like ConA's subunit-subunit interfaces): the
    # correction should be small.
    s_cav = sasa_particles(prot; exclude_cavities=true)
    @test sasa(s_cav) <= sasa(s_default) # can only remove area, never add
    @test sasa(s_cav) ≈ sasa(s_default) rtol = 0.01

    # Results should be essentially insensitive to the exact cutoff multiple over a
    # wide, physically sensible range (this is the point of connecting on the dot
    # cloud itself rather than on a resampled grid).
    s_06 = sasa(sasa_particles(prot; exclude_cavities=true, cavity_dot_cutoff=0.6))
    s_13 = sasa(sasa_particles(prot; exclude_cavities=true, cavity_dot_cutoff=1.3))
    @test s_06 ≈ s_13 rtol = 1e-3

    # unitcell (periodic boundary conditions) is not supported together with cavity
    # exclusion; must fail loudly rather than silently ignore one or the other.
    @test_throws "not currently supported together with" sasa_particles(
        prot; exclude_cavities=true, unitcell=[100.0, 100.0, 100.0],
    )
end

@testitem "exclude_cavities: reproduces SurfaceRacer on 3CNA" begin
    using PDBTools

    # Regression pins from running SurfaceRacer 5.0 itself (Tsodikov, Record & Sergeev,
    # 2002), locally, with its "1 - Richards (1977)" radii option and a 1.4 Å probe, in
    # ASA-only mode, on the exact same biological tetramer assembly and dimer split used
    # in the "RichardsUnitedAtomRadii validated against SurfaceRacer" testitem
    # (richards.jl): tetramer = 34052.06 Å², each dimer = 19469.37 Å² (ΔASA = 4886.7 Å²,
    # matching Table S5 of Knowles et al. 2015 to within rounding).
    #
    # The dot-graph cavity exclusion implemented in this file gets much closer to that
    # than the voxel-grid approach it replaced (see this file's header comment): the
    # *difference* (2*dimer - tetramer) -- the quantity that actually enters an m-value
    # -- lands within ~1% of SurfaceRacer's ΔASA, essentially independent of the exact
    # cutoff chosen (tested 0.6-1.3 Å; all agree to <0.01%). This is not an exact,
    # bit-for-bit reproduction of SurfaceRacer's analytical algorithm (individual
    # tetramer/dimer totals still sit ~1-3% low), but it is a large, robust improvement
    # over the uncorrected default (4318.8 Å², ~12% short of target).
    cna = wget("3CNA", "protein"; assembly=1)
    d1 = select(cna, at -> chain(at) in ("A", "A-2"))
    d2 = select(cna, at -> chain(at) in ("A-3", "A-4"))

    cavkw = (exclude_cavities=true, cavity_dot_cutoff=1.0)
    s_tet = sasa(sasa_particles(RichardsUnitedAtomRadii, cna; radii_set=:set1, cavkw...))
    s_d1 = sasa(sasa_particles(RichardsUnitedAtomRadii, d1; radii_set=:set1, cavkw...))
    s_d2 = sasa(sasa_particles(RichardsUnitedAtomRadii, d2; radii_set=:set1, cavkw...))

    @test s_tet ≈ 33771.6 rtol = 1e-3
    @test s_d1 ≈ 19329.8 rtol = 1e-3
    @test s_d2 ≈ 19283.0 rtol = 1e-3
    delta = s_d1 + s_d2 - s_tet
    @test delta ≈ 4841.2 rtol = 1e-3
    @test delta ≈ 4886.7 rtol = 0.02   # within ~1% of SurfaceRacer's own ΔASA
    @test delta > 4318.8               # a real improvement over the uncorrected default

    # Robustness: the *difference* should barely move across a wide cutoff range, even
    # though it necessarily requires re-detecting connectivity from scratch each time.
    for cutoff in (0.6, 0.8, 1.3)
        cavkw2 = (exclude_cavities=true, cavity_dot_cutoff=cutoff)
        d1v = sasa(sasa_particles(RichardsUnitedAtomRadii, d1; radii_set=:set1, cavkw2...))
        d2v = sasa(sasa_particles(RichardsUnitedAtomRadii, d2; radii_set=:set1, cavkw2...))
        tetv = sasa(sasa_particles(RichardsUnitedAtomRadii, cna; radii_set=:set1, cavkw2...))
        @test (d1v + d2v - tetv) ≈ delta rtol = 1e-3
    end
end
