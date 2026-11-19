export helical_wheel, HelicalWheel, hydrophobic_moment

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

"""
    hydrophobic_moment(sequence; hydrophobicity = Fauchère-Pliska scale)

Compute the mean hydrophobicity and the hydrophobic moment (Eisenberg) of an
α-helix with the given `sequence` (one-letter codes), assuming 100° between
consecutive residues.

The default hydrophobicity scale is that of Fauchère & Pliska (1983). Another scale
can be given as a dictionary from one-letter codes (`Char`) to values.

Returns a named tuple `(mean_hydrophobicity, μ, μx, μy)`, where `μ = √(μx² + μy²)`
is the moment per residue.

# Example

```julia-repl
julia> using PDBTools

julia> m = hydrophobic_moment("LKKLLKLLKKLLKL");

julia> round(m.μ; digits = 2)
0.81
```
"""
function hydrophobic_moment(sequence::AbstractString; hydrophobicity::AbstractDict = FAUCHERE_PLISKA)
    seq = collect(uppercase(sequence))
    n = length(seq)
    n > 0 || throw(ArgumentError("The sequence is empty."))
    unknown = unique(filter(aa -> !haskey(hydrophobicity, aa), seq))
    isempty(unknown) ||
        throw(ArgumentError("Residues without a value in the hydrophobicity scale: $(join(unknown, ", "))"))
    angles = _helical_wheel_angles(n)
    h = [hydrophobicity[aa] for aa in seq]
    μx = sum(h[i] * cos(angles[i]) for i in 1:n) / n
    μy = sum(h[i] * sin(angles[i]) for i in 1:n) / n
    return (mean_hydrophobicity = sum(h) / n, μ = sqrt(μx^2 + μy^2), μx = μx, μy = μy)
end

"""
    HelicalWheel

Helical wheel of a helix, with one value per residue, as returned by [`helical_wheel`](@ref).
Plot it with `plot(wheel)` (requires `using Plots`).

# Fields
- `sequence`: sequence of the helix (one-letter codes).
- `values`: one value per residue.
- `residue_numbers`: the residue numbers.
- `mean_hydrophobicity`, `μ`, `μx`, `μy`: mean hydrophobicity and hydrophobic moment of
  the helix (see [`hydrophobic_moment`](@ref)).
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
`hydrophobicity` (default: Fauchère & Pliska, 1983; see [`hydrophobic_moment`](@ref)).

# Examples

```julia-repl
julia> using PDBTools, LazyArtifacts, Plots

julia> wheel = helical_wheel("LKKLLKLLKKLLKL", rand(14));

julia> plot(wheel; colorbar_title = "property")

julia> pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> helix = select(pdb, "resnum >= 10 and resnum <= 25");

julia> plot(helical_wheel(helix; aggregate = sum); colorbar_title = "CN")
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
        moment.mean_hydrophobicity, moment.μ, moment.μx, moment.μy,
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
    m = hydrophobic_moment(repeat("L", 18))
    @test m.μ ≈ 0.0 atol = 1e-12
    @test m.mean_hydrophobicity ≈ PDBTools.FAUCHERE_PLISKA['L']
    # an amphipathic sequence has a non-zero moment; lowercase is accepted
    @test hydrophobic_moment("lkkllkllkkllkl").μ ≈ 0.8132 atol = 1e-3
    @test hydrophobic_moment("LKKLLKLLKKLLKL").μ ≈ hydrophobic_moment("lkkllkllkkllkl").μ
    # custom scale
    @test hydrophobic_moment("LK"; hydrophobicity = Dict('L' => 1.0, 'K' => 0.0)).mean_hydrophobicity ≈ 0.5
    @test_throws ArgumentError hydrophobic_moment("LKXL")
    @test_throws ArgumentError hydrophobic_moment("")
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
    @test (w.mean_hydrophobicity, w.μ, w.μx, w.μy) == (m.mean_hydrophobicity, m.μ, m.μx, m.μy)
    @test parse_show(w) ≈ """
        HelicalWheel of 14 residues (11-24): LKKLLKLLKKLLKL
            Mean hydrophobicity: $(m.mean_hydrophobicity)
            Hydrophobic moment: $(m.μ)
        """
    @test_throws ArgumentError helical_wheel(seq, ones(n - 1))
    @test_throws ArgumentError helical_wheel(seq, ones(n); residue_numbers = 1:n-1)
    @test_throws ArgumentError helical_wheel("LKXL", ones(4))

    # from a selection of atoms: the values are aggregate(property, residue)
    pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    atoms = select(pdb, at -> resnum(at) <= 18)
    residues = collect(eachresidue(atoms))
    w = helical_wheel(atoms; property = mass, aggregate = sum)
    @test w.values ≈ [sum(mass, r) for r in residues]
    @test w.sequence == join(oneletter.(resname.(residues)))
    @test w.residue_numbers == resnum.(residues)
    w = helical_wheel(atoms; property = mass)
    @test w.values ≈ [sum(mass, r) / length(r) for r in residues]
    @test_throws ArgumentError helical_wheel(atoms[1:0])
    @test_throws ArgumentError helical_wheel(read_pdb(PDBTools.test_file("dimer.pdb"), "protein"))
end
