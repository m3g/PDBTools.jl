export helical_wheel, HelicalWheel, hydrophobic_moment, mean_hydrophobicity

# Hydrophobicity scale of Fauchère & Pliska (1983), indexed by the one-letter residue code.
# Default scale of hydrophobic_moment and helical_wheel.
const FAUCHERE_PLISKA = Dict(
    'A' => 0.31,  'R' => -1.01, 'N' => -0.60, 'D' => -0.77,
    'C' => 1.54,  'Q' => -0.22, 'E' => -0.64, 'G' => 0.00,
    'H' => 0.13,  'I' => 1.80,  'L' => 1.70,  'K' => -0.99,
    'M' => 1.23,  'F' => 1.79,  'P' => 0.72,  'S' => -0.04,
    'T' => 0.26,  'W' => 2.25,  'Y' => 0.96,  'V' => 1.22,
)

# Angle between consecutive residues in an α-helix, in degrees
const HELIX_ANGLE_STEP = 100

# Angular positions (radians) of the residues in the helical wheel
_helical_wheel_angles(n::Integer) = [deg2rad(-i * HELIX_ANGLE_STEP % 360) for i in 0:n-1]

# Mean of f over the elements of x (Statistics is not a dependency)
_mean(f, x) = sum(f, x) / length(x)

# Hydrophobicity of each residue of the sequence, in the given scale
function _hydrophobicities(sequence::AbstractString, hydrophobicity::AbstractDict)
    seq = collect(uppercase(sequence))
    isempty(seq) && throw(ArgumentError("The sequence is empty."))
    unknown = unique(filter(aa -> !haskey(hydrophobicity, aa), seq))
    isempty(unknown) ||
        throw(ArgumentError("Residues without a value in the hydrophobicity scale: $(join(unknown, ", "))"))
    return [Float64(hydrophobicity[aa]) for aa in seq]
end

"""
    hydrophobic_moment(sequence; hydrophobicity = Fauchère-Pliska scale)

Compute the hydrophobic moment (Eisenberg) of an α-helix with the given `sequence`
(one-letter codes), assuming 100° between consecutive residues.

Returns the moment per residue, as a 2D vector `(μx, μy)`, in the plane perpendicular
to the axis of the helix. Its magnitude, `hypot(moment...)` (or `norm(moment)`, from
`LinearAlgebra`), is the amphipathicity of the helix.

The default hydrophobicity scale is that of Fauchère & Pliska (1983). Another scale
can be given as a dictionary from one-letter codes (`Char`) to values.

See also [`mean_hydrophobicity`](@ref).

# Example

```julia-repl
julia> using PDBTools

julia> m = hydrophobic_moment("LKKLLKLLKKLLKL");

julia> round(hypot(m...); digits = 2)
0.81
```
"""
function hydrophobic_moment(sequence::AbstractString; hydrophobicity::AbstractDict = FAUCHERE_PLISKA)
    h = _hydrophobicities(sequence, hydrophobicity)
    n = length(h)
    angles = _helical_wheel_angles(n)
    μx = sum(h[i] * cos(angles[i]) for i in 1:n) / n
    μy = sum(h[i] * sin(angles[i]) for i in 1:n) / n
    return SVector(μx, μy)
end

"""
    mean_hydrophobicity(sequence; hydrophobicity = Fauchère-Pliska scale)

Compute the mean hydrophobicity of the residues of the given `sequence` (one-letter codes).

The default hydrophobicity scale is that of Fauchère & Pliska (1983). Another scale
can be given as a dictionary from one-letter codes (`Char`) to values.

See also [`hydrophobic_moment`](@ref).

# Example

```julia-repl
julia> using PDBTools

julia> mean_hydrophobicity("LKKLLKLLKKLLKL")
0.5471428571428572
```
"""
mean_hydrophobicity(sequence::AbstractString; hydrophobicity::AbstractDict = FAUCHERE_PLISKA) =
    _mean(identity, _hydrophobicities(sequence, hydrophobicity))

"""
    HelicalWheel

Helical wheel of a helix, with one value per residue, as returned by [`helical_wheel`](@ref).
Plot it with `plot(wheel)` (requires `using Plots`).

# Fields
- `sequence`: sequence of the helix (one-letter codes).
- `values`: one value per residue.
- `residue_numbers`: the residue numbers.
- `mean_hydrophobicity`: the mean hydrophobicity of the helix (see [`mean_hydrophobicity`](@ref)).
- `μ`, `μx`, `μy`: the magnitude and the components of the hydrophobic moment of the helix
  (see [`hydrophobic_moment`](@ref)).
"""
struct HelicalWheel
    sequence::String
    values::Vector{Float64}
    residue_numbers::Vector{Int}
    mean_hydrophobicity::Float64
    μ::Float64
    μx::Float64
    μy::Float64
end

function Base.show(io::IO, ::MIME"text/plain", w::HelicalWheel)
    print(io, chomp("""
    HelicalWheel of $(length(w.sequence)) residues ($(first(w.residue_numbers))-$(last(w.residue_numbers))): $(w.sequence)
        Mean hydrophobicity: $(w.mean_hydrophobicity)
        Hydrophobic moment: $(w.μ)
    """))
end

"""
    helical_wheel(sequence, values; residue_numbers = 1:length(sequence), hydrophobicity)
    helical_wheel(atoms; property = beta, aggregate = mean, hydrophobicity)

Compute the helical wheel of a helix, with one value per residue. The values can be any
per-residue property (coordination number, SASA, RMSF, interaction energy, conservation, etc.),
in the same order as the sequence.

Returns a [`HelicalWheel`](@ref), which can be plotted with `plot(wheel)` (requires `using Plots`).

With a vector of atoms (of a single chain), the sequence, the residue numbers and the values
are taken from the structure: the value of each residue is `aggregate(property, residue)`,
that is, `property` is computed for each atom and combined for the residue with `aggregate`:

- `property`: function applied to each atom (default: `beta`; can be `occup`, `charge`, etc.).
- `aggregate`: `sum` for additive properties (e.g. atomic contributions to a coordination
  number), or the mean (default) for intensive properties (B-factor, RMSF). Any function with
  the `aggregate(f, collection)` form can be used, as `sum`, `maximum`, `minimum`, or `mean`
  (from `Statistics`).

The mean hydrophobicity and the hydrophobic moment are computed with the scale given by
`hydrophobicity` (default: Fauchère & Pliska, 1983; see [`mean_hydrophobicity`](@ref) and
[`hydrophobic_moment`](@ref)).

# Examples

```julia-repl
julia> using PDBTools, LazyArtifacts, Plots

julia> wheel = helical_wheel("LKKLLKLLKKLLKL", rand(14));

julia> plot(wheel; colorbar_title = "property")

julia> pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> helix = select(pdb, "resnum >= 13 and resnum <= 29"); # an α-helix

julia> plot(helical_wheel(helix; property = mass, aggregate = sum); colorbar_title = "mass")
```

Values from other packages can be passed directly. For example, with a
`ResidueContributions` object `rc` of ComplexMixtures.jl, the contribution of
each residue at distance `d` is:

```julia
id = findfirst(>=(d), rc.d)
values = [rc.residue_contributions[i][id] for i in eachindex(rc.residue_contributions)]
wheel = helical_wheel(sequence, values)
```
"""
function helical_wheel(
    sequence::AbstractString,
    values::AbstractVector{<:Real};
    residue_numbers::AbstractVector{<:Integer} = 1:length(sequence),
    hydrophobicity::AbstractDict = FAUCHERE_PLISKA,
)
    n = length(sequence)
    length(values) == n ||
        throw(ArgumentError("values has $(length(values)) elements, but the sequence has $n residues."))
    length(residue_numbers) == n ||
        throw(ArgumentError("residue_numbers has $(length(residue_numbers)) elements, but the sequence has $n residues."))
    # Mean hydrophobicity and hydrophobic moment (also validates the sequence)
    moment = hydrophobic_moment(sequence; hydrophobicity)
    return HelicalWheel(
        uppercase(sequence), float.(values), collect(residue_numbers),
        mean_hydrophobicity(sequence; hydrophobicity), norm(moment), moment[1], moment[2],
    )
end

function helical_wheel(
    atoms::AbstractVector{<:Atom};
    property = beta,
    aggregate = _mean,
    hydrophobicity::AbstractDict = FAUCHERE_PLISKA,
)
    isempty(atoms) && throw(ArgumentError("The atom selection is empty."))
    chains = unique(chain.(atoms))
    length(chains) == 1 ||
        throw(ArgumentError("The selection contains more than one chain ($(join(chains, ", "))). Select a single chain."))
    residues = collect(eachresidue(atoms))
    values = [aggregate(property, r) for r in residues]
    sequence = join(oneletter.(resname.(residues)))
    return helical_wheel(sequence, values; residue_numbers = resnum.(residues), hydrophobicity)
end

@testitem "hydrophobic_moment" begin
    using PDBTools
    # 18 residues at 100° complete exactly 5 turns: a uniform sequence has zero moment
    using StaticArrays: SVector
    norm(m) = hypot(m...)
    m = hydrophobic_moment(repeat("L", 18))
    @test m isa SVector{2,Float64}
    @test norm(m) ≈ 0.0 atol = 1e-12
    @test mean_hydrophobicity(repeat("L", 18)) ≈ PDBTools.FAUCHERE_PLISKA['L']
    # an amphipathic sequence has a non-zero moment; lowercase is accepted
    @test norm(hydrophobic_moment("lkkllkllkkllkl")) ≈ 0.8132 atol = 1e-3
    @test hydrophobic_moment("LKKLLKLLKKLLKL") ≈ hydrophobic_moment("lkkllkllkkllkl")
    @test mean_hydrophobicity("LKKLLKLLKKLLKL") ≈ (8 * 1.70 + 6 * -0.99) / 14
    # custom scale
    scale = Dict('L' => 1.0, 'K' => 0.0)
    @test mean_hydrophobicity("LK"; hydrophobicity = scale) ≈ 0.5
    @test norm(hydrophobic_moment("LLLLLLLLLLLLLLLLLL"; hydrophobicity = scale)) ≈ 0.0 atol = 1e-12
    for f in (hydrophobic_moment, mean_hydrophobicity)
        @test_throws ArgumentError f("LKXL")
        @test_throws ArgumentError f("")
    end
end

@testitem "helical_wheel" begin
    using PDBTools, LazyArtifacts
    using ShowMethodTesting
    seq = "LKKLLKLLKKLLKL"
    n = length(seq)
    w = helical_wheel(lowercase(seq), collect(1.0:n); residue_numbers = 11:24)
    @test w isa HelicalWheel
    @test w.sequence == seq
    @test w.values == collect(1.0:n)
    @test w.residue_numbers == collect(11:24)
    m = hydrophobic_moment(seq)
    @test w.mean_hydrophobicity == mean_hydrophobicity(seq)
    @test (w.μx, w.μy) == (m[1], m[2])
    @test w.μ ≈ sqrt(m[1]^2 + m[2]^2)
    @test parse_show(w) ≈ """
        HelicalWheel of 14 residues (11-24): LKKLLKLLKKLLKL
            Mean hydrophobicity: $(w.mean_hydrophobicity)
            Hydrophobic moment: $(w.μ)
        """
    @test_throws ArgumentError helical_wheel(seq, ones(n - 1))
    @test_throws ArgumentError helical_wheel(seq, ones(n); residue_numbers = 1:n-1)
    @test_throws ArgumentError helical_wheel("LKXL", ones(4))

    # from a selection of atoms: the values are aggregate(property, residue)
    pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    atoms = select(pdb, at -> resnum(at) <= 18)
    residues = collect(eachresidue(atoms))
    w = helical_wheel(atoms) # default: mean B-factor of each residue
    @test w.values ≈ [sum(beta, r) / length(r) for r in residues]
    @test !allequal(w.values)
    @test w.sequence == join(oneletter.(resname.(residues)))
    @test w.residue_numbers == resnum.(residues)
    @test helical_wheel(atoms; aggregate = sum).values ≈ [sum(beta, r) for r in residues]
    @test helical_wheel(atoms; property = mass, aggregate = sum).values ≈ [sum(mass, r) for r in residues]
    @test_throws ArgumentError helical_wheel(atoms[1:0])
    @test_throws ArgumentError helical_wheel(read_pdb(PDBTools.test_file("dimer.pdb"), "protein"))
end
