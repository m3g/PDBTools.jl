export helical_wheel, residue_values, hydrophobic_moment

"""
    FAUCHERE_PLISKA

Hydrophobicity scale of Fauchère & Pliska (1983), as a `Dict{Char, Float64}`
indexed by the one-letter residue code. Default scale of [`hydrophobic_moment`](@ref)
and [`helical_wheel`](@ref).
"""
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

_mean(x) = sum(x) / length(x)

"""
    hydrophobic_moment(sequence; hydrophobicity = PDBTools.FAUCHERE_PLISKA)

Compute the mean hydrophobicity and the hydrophobic moment (Eisenberg) of an
α-helix with the given `sequence` (one-letter codes), assuming 100° between
consecutive residues.

Returns a named tuple `(mean_hydro, μ, μx, μy)`, where `μ = √(μx² + μy²)`
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
    return (mean_hydro = sum(h) / n, μ = sqrt(μx^2 + μy^2), μx = μx, μy = μy)
end

"""
    residue_values(atoms; aggregate = mean, property = beta)

Group the atoms by residue and reduce `property` of the atoms of each residue
with `aggregate`. Returns `(values, sequence, residue_numbers)`.

- `aggregate`: any function that takes a vector, e.g. `sum`, `maximum`, `minimum`.
  The default is the mean. Use `sum` for additive quantities (e.g. atomic
  contributions to the coordination number) and the mean for intensive
  properties (B-factor, RMSF).
- `property`: function applied to each atom (default: `beta`; can be `occup`,
  `charge`, etc.).

The selection must contain a single chain.

# Example

```julia-repl
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein and resnum <= 10");

julia> values, sequence, numbers = residue_values(atoms; aggregate = sum);
```
"""
function residue_values(atoms::AbstractVector{<:Atom}; aggregate = _mean, property = beta)
    isempty(atoms) && throw(ArgumentError("The atom selection is empty."))
    chains = unique(chain.(atoms))
    length(chains) == 1 ||
        throw(ArgumentError("The selection contains more than one chain ($(join(chains, ", "))). Select a single chain."))
    residues = collect(eachresidue(atoms))
    values = [aggregate(property.(res)) for res in residues]
    sequence = join(oneletter.(resname.(residues)))
    numbers = resnum.(residues)
    return values, sequence, numbers
end

"""
    helical_wheel(sequence, values; kwargs...)
    helical_wheel(atoms; aggregate = mean, property = beta, kwargs...)

Draw a helical wheel with each residue colored by `values[i]`.
`values` can be any per-residue property (coordination number, SASA, RMSF,
interaction energy, conservation, etc.), in the same order as the sequence.

!!! note
    This function requires `using Plots`.

With a vector of atoms, the sequence, the residue numbers and the values are
taken from the structure with [`residue_values`](@ref).

# Arguments
- `sequence::AbstractString`: sequence of the helix (one-letter codes).
- `values::AbstractVector{<:Real}`: one value per residue.

# Keyword arguments
- `residue_numbers = 1:length(sequence)`: actual residue numbers (used only in
  the labels of the first and last residues).
- `normalize = true`: apply min-max normalization to `[0, 1]`.
- `colormap = :tempo`: color map of the values.
- `colorbar_title = ""`: title of the color bar.
- `clims = nothing`: fixed limits of the color scale, e.g. `(0, 1)`.
  Use to compare several helices on the same scale.
- `hydrophobicity = PDBTools.FAUCHERE_PLISKA`: scale used in the hydrophobic moment.
- `show_moment = true`: draw the arrow of the hydrophobic moment.

Returns `(plot, (mean_hydro = ..., μ = ...))`.

# Examples

```julia-repl
julia> using PDBTools, Plots

julia> pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> helix = select(pdb, "chain A and resnum >= 10 and resnum <= 25");

julia> p, info = helical_wheel(helix; aggregate = sum, colorbar_title = "CN");
```

Values from other packages can be passed directly. For example, with a
`ResidueContributions` object `rc` of ComplexMixtures.jl, the contribution of
each residue at distance `d` is:

```julia
id = findfirst(>=(d), rc.d)
values = [rc.residue_contributions[i][id] for i in eachindex(rc.residue_contributions)]
p, info = helical_wheel(sequence, values)
```
"""
function helical_wheel end

@testitem "hydrophobic_moment" begin
    using PDBTools
    # 18 residues at 100° complete exactly 5 turns: a uniform sequence has zero moment
    m = hydrophobic_moment(repeat("L", 18))
    @test m.μ ≈ 0.0 atol = 1e-12
    @test m.mean_hydro ≈ PDBTools.FAUCHERE_PLISKA['L']
    # an amphipathic sequence has a non-zero moment; lowercase is accepted
    @test hydrophobic_moment("lkkllkllkkllkl").μ ≈ 0.8132 atol = 1e-3
    @test hydrophobic_moment("LKKLLKLLKKLLKL").μ ≈ hydrophobic_moment("lkkllkllkkllkl").μ
    @test_throws ArgumentError hydrophobic_moment("LKXL")
    @test_throws ArgumentError hydrophobic_moment("")
end

@testitem "residue_values" begin
    using PDBTools, LazyArtifacts
    pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    ch = chain(first(pdb))
    atoms = select(pdb, at -> chain(at) == ch && resnum(at) <= resnum(first(pdb)) + 9)
    values, sequence, numbers = residue_values(atoms)
    residues = collect(eachresidue(atoms))
    @test length(values) == length(sequence) == length(numbers) == length(residues)
    @test values[1] ≈ sum(beta.(residues[1])) / length(residues[1])
    values_sum, _, _ = residue_values(atoms; aggregate = sum)
    @test values_sum[1] ≈ sum(beta.(residues[1]))
    @test_throws ArgumentError residue_values(atoms[1:0])
end
