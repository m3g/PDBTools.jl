export secondary_structure_diagram, SecondaryStructureDiagram

"""
    SecondaryStructureDiagram

Secondary structure of a protein, with (optionally) one value per residue, as computed by
[`secondary_structure_diagram`](@ref). It is plotted with `plot(diagram)` (requires a Makie backend, e.g.
`using CairoMakie`).

# Fields
- `ss::Vector{SSData}`: secondary structure of each residue (see [`dssp_run`](@ref)).
- `values::Union{Nothing,Vector{Float64}}`: value of each residue, or `nothing`.

"""
struct SecondaryStructureDiagram
    ss::Vector{SSData}
    values::Union{Nothing,Vector{Float64}}
end

function Base.show(io::IO, ::MIME"text/plain", d::SecondaryStructureDiagram)
    chains = unique(s.chain for s in d.ss)
    print(io, "SecondaryStructureDiagram of $(length(d.ss)) residues ($(length(chains)) chain$(length(chains) > 1 ? "s" : "")), ")
    if isnothing(d.values)
        print(io, "without values")
    else
        lo, hi = round.(extrema(d.values); sigdigits = 4)
        print(io, "with values in [$lo, $hi]")
    end
end

"""
    secondary_structure_diagram(atoms; property = beta, aggregate = mean, ss_method = :dssp)
    secondary_structure_diagram(atoms, values; ss_method = :dssp)
    secondary_structure_diagram(ss, values = nothing)

Compute the secondary structure diagram of a protein: the secondary structure of each
residue, with one value per residue. Returns an [`SecondaryStructureDiagram`](@ref), which is plotted
with `plot(diagram)` (requires a Makie backend, e.g. `using CairoMakie`).

DSSP and STRIDE need the backbone atoms of each residue, so coarse-grained structures
(e.g. Martini or SIRAH) are not supported: all residues are assigned as loops.

# Input

- `atoms`: the secondary structure is computed with [`dssp_run`](@ref) (`ss_method = :dssp`,
  default) or [`stride_run`](@ref) (`ss_method = :stride`). The value of each residue is
  `aggregate(property, residue)`: `property` is computed for each atom, and combined for the
  residue with `aggregate`:
  - `property`: function applied to each atom (default: `beta`; can be `occup`, `charge`,
    or any function of the atom, e.g. `at -> beta(at)^2`).
  - `aggregate`: the mean (default) for intensive properties (B-factors, RMSF), or `sum` for
    additive properties (e.g. atomic contributions to a coordination number). Any function
    with the `aggregate(f, collection)` form can be used, as `maximum` or `minimum`.
- `atoms, values`: one value per residue, in the order of `eachresidue(atoms)`.
- `ss`: the output of [`dssp_run`](@ref) or [`stride_run`](@ref), such that the secondary
  structure is not computed again, and optionally the `values`, one per element of `ss`.
  Without values, only the secondary structure is drawn.

Only protein residues have a secondary structure: the values of other residues are ignored.

# Examples

```julia-repl
julia> using PDBTools, LazyArtifacts, CairoMakie

julia> atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> d = secondary_structure_diagram(atoms) # colored by the mean B-factor of each residue
SecondaryStructureDiagram of 104 residues (1 chain), with values in [13.1, 26.36]

julia> plot(d; colorbar_title = "B-factor")

julia> ss = dssp_run(atoms);

julia> plot(secondary_structure_diagram(ss)) # only the secondary structure
```

"""
function secondary_structure_diagram(
    atoms::AbstractVector{<:Atom};
    property = beta,
    aggregate = _mean,
    ss_method::Symbol = :dssp,
)
    residues = collect(eachresidue(atoms))
    isempty(residues) && throw(ArgumentError("The atom selection has no residues."))
    property isa Function || throw(ArgumentError(
        "property must be a function of the atom (e.g. beta or occup), got $(typeof(property))."))
    # If the property is constant within each residue, it is probably already a
    # per-residue value; summing would multiply it by the number of atoms.
    if aggregate === sum
        per_atom = [property.(r) for r in residues]
        # Base.all: PDBTools defines its own `all` (the "all" selection)
        if Base.all(v -> Base.all(≈(first(v)), v), per_atom) && any(v -> length(v) > 1, per_atom)
            @warn "aggregate = sum, but the property is constant within each residue " *
                  "(it looks like a per-residue value). Consider the mean."
        end
    end
    values = [Float64(aggregate(property, r)) for r in residues]
    return _secondary_structure_diagram(atoms, residues, values, ss_method)
end

function secondary_structure_diagram(atoms::AbstractVector{<:Atom}, values::AbstractVector{<:Real}; ss_method::Symbol = :dssp)
    residues = collect(eachresidue(atoms))
    isempty(residues) && throw(ArgumentError("The atom selection has no residues."))
    length(values) == length(residues) || throw(ArgumentError(
        "values has $(length(values)) elements, but there are $(length(residues)) residues."))
    return _secondary_structure_diagram(atoms, residues, Float64.(values), ss_method)
end

function _secondary_structure_diagram(atoms, residues, values, ss_method)
    ss = if ss_method == :dssp
        dssp_run(atoms)
    elseif ss_method == :stride
        stride_run(atoms)
    else
        throw(ArgumentError("ss_method must be :dssp or :stride, got :$ss_method"))
    end
    # The output of dssp_run and stride_run has one element per protein residue, in order
    # (only of the first model, for DSSP): non-protein residues have no secondary structure.
    ss_values = values[findall(r -> isprotein(first(r)), residues)][1:length(ss)]
    return SecondaryStructureDiagram(ss, ss_values)
end

function secondary_structure_diagram(ss::AbstractVector{<:SSData}, values::Union{Nothing,AbstractVector{<:Real}} = nothing)
    isempty(ss) && throw(ArgumentError("The secondary structure data is empty."))
    if !isnothing(values)
        length(values) == length(ss) || throw(ArgumentError(
            "values has $(length(values)) elements, but there are $(length(ss)) residues."))
    end
    return SecondaryStructureDiagram(collect(ss), isnothing(values) ? nothing : Float64.(values))
end

@testitem "secondary_structure_diagram" begin
    using PDBTools, LazyArtifacts
    atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    residues = collect(eachresidue(atoms))
    nres = length(residues)

    d = secondary_structure_diagram(atoms)
    @test d isa SecondaryStructureDiagram
    @test length(d.ss) == nres
    @test d.values ≈ [sum(beta, r) / length(r) for r in residues]
    @test !allequal(d.values) # the test file has non-trivial B-factors
    @test secondary_structure_diagram(atoms; aggregate = sum).values ≈ [sum(beta, r) for r in residues]
    @test secondary_structure_diagram(atoms; property = occup, aggregate = maximum).values ≈ [maximum(occup, r) for r in residues]
    @test secondary_structure_diagram(atoms; property = at -> 2 * beta(at), aggregate = sum).values ≈ [2 * sum(beta, r) for r in residues]
    @test length(secondary_structure_diagram(atoms; ss_method = :stride).ss) == nres
    values = collect(range(-1.0, 1.0; length = nres))
    @test secondary_structure_diagram(atoms, values).values == values
    @test occursin("SecondaryStructureDiagram of 104 residues (1 chain), with values", sprint(show, MIME"text/plain"(), d))

    # from the output of dssp_run
    ss = dssp_run(atoms)
    @test isnothing(secondary_structure_diagram(ss).values)
    @test occursin("without values", sprint(show, MIME"text/plain"(), secondary_structure_diagram(ss)))
    @test secondary_structure_diagram(ss, 1:length(ss)).values == collect(1.0:length(ss))

    # non-protein residues are ignored
    pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein or resname SOD")
    @test length(secondary_structure_diagram(pdb).ss) == nres

    # warning when aggregate = sum but the property is constant within each residue
    atoms_const = deepcopy(atoms)
    for (i, r) in enumerate(eachresidue(atoms_const)), at in r
        at.beta = i
    end
    @test_logs (:warn, r"aggregate = sum") match_mode = :any secondary_structure_diagram(atoms_const; aggregate = sum)

    @test_throws ArgumentError secondary_structure_diagram(atoms; ss_method = :foo)
    @test_throws ArgumentError secondary_structure_diagram(atoms; property = :beta)
    @test_throws ArgumentError secondary_structure_diagram(atoms, ones(nres - 1))
    @test_throws ArgumentError secondary_structure_diagram(atoms[1:0])
    @test_throws ArgumentError secondary_structure_diagram(ss[1:0])
    @test_throws ArgumentError secondary_structure_diagram(ss, ones(length(ss) - 1))
end
