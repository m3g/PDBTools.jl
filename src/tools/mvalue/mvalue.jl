export mvalue
export MoeserHorinek, AutonBolen

abstract type MvalueModel end
struct MoeserHorinek <: MvalueModel end
struct AutonBolen <: MvalueModel end

include("./data.jl")

struct MValue{T}
    nresidues::Int
    tot::Float32
    bb::Float32
    sc::Float32
    residue_contributions_bb::Vector{Float32}
    residue_contributions_sc::Vector{Float32}
end
function Base.show(io::IO, ::MIME"text/plain", m::MValue)
    print(io, chomp("""
    $(typeof(m)) - $(m.nresidues) residues.
        Total m-value: $(m.tot) kcal mol⁻¹
        Backbone contributions: $(m.bb) kcal mol⁻¹
        Side-chain contributions: $(m.sc) kcal mol⁻¹
    """))
end

"""
    mvalue(
        sasa_initial::SASA{3,<:AbstractVector{<:Atom}},    
        sasa_final::SASA{3,<:AbstractVector{<:Atom}},
        cosolvent::String;
        sel::Union{String,Function}=all,
        model::Type{<:MvalueModel}=AutonBolen,
        backbone::Function = isbackbone,
        sidechain::Function = issidechain,
    )

Calculates the m-value (transfer free energy of a protein in 1M solution, in `kcal/mol`) using the Tanford transfer model,
as implemented by Moeser and Horinek [1] or by Auton and Bolen [2,3].

# Positional Arguments

- `sasa_initial::SASA{3,<:AbstractVector{<:Atom}}`: SASA object representing the initial state (e.g., native state).
- `sasa_final::SASA{3,<:AbstractVector{<:Atom}}`: SASA object representing the final state (e.g., denatured state).
- `cosolvent::String`: The cosolvent used in the solution. One of $(join('"' .* sort!(unique(keys(PDBTools.cosolvent_column)) .* '"'; by=lowercase),", ")) 

# Keyword Arguments

- `sel::Union{String,Function}=all`: Selection of atoms to consider in the calculation. Can be a selection string or a function that takes an `Atom` and returns a `Bool`.
- `model::Type{<:MvalueModel}=AutonBolen`: The model to use for the calculation. Either `MoeserHorinek` or `AutonBolen`.
- `backbone::Function = PDBTools.isbackbone`: Function to identify backbone atoms.
- `sidechain::Function = PDBTools.issidechain`: Function to identify side chain atoms.

# Returns

A `MValue` object, with fields:

- `ntatoms::Int`: Number of atoms considered.
- `tot::Float32`: Total m-value (kcal/mol/M).
- `bb::Float32`: Backbone contribution to the m-value (kcal/mol/M).
- `sc::Float32`: Side chain contribution to the m-value (kcal/mol/M).
- `residue_contributions_bb::Vector{Float32}`: Backbone contributions of each residue to the m-value.
- `residue_contributions_sc::Vector{Float32}`: Side-chain contributions of each residue to the m-value.

## Example

```julia
using PDBTools

initial_state = read_pdb("native.pdb")
final_state = read_pdb("desnat.pdb")

sasa_initial = sasa_particles(initial_state)
sasa_final = sasa_particles(final_state)

mvalue(sasa_initial, sasa_final, "chain A"; model=AutonBolen, cosolvent="TMAO")
```

## References

1. https://doi.org/10.1021/acs.jpcb.7b02138
2. https://doi.org/10.1016/s0076-6879(07)28023-1
3. https://www.pnas.org/doi/10.1073/pnas.0706251104

"""
function mvalue(
    sasa_initial::SASA{3,<:AbstractVector{<:Atom}},
    sasa_final::SASA{3,<:AbstractVector{<:Atom}},
    cosolvent::String;
    sel::Union{String,Function}=all,
    model::Type{<:MvalueModel}=AutonBolen,
    backbone::F1=isbackbone,
    sidechain::F2=issidechain,
) where {F1<:Function, F2<:Function}
    selector = Select(sel)
    ats_initial = filter(selector, sasa_initial.particles)
    ats_final = filter(selector, sasa_final.particles)
    residues_initial = collect(eachresidue(ats_initial))
    residues_final = collect(eachresidue(ats_final))
    if length(residues_initial) != length(residues_final)
        throw(ArgumentError("""\n
            Initial and final states do not have the same number of residues.
            Got $(length(residues_initial)) and $(length(residues_final)) residues.
        
        """
        ))
    end
    residue_contributions_bb = zeros(Float32, length(residues_initial))
    residue_contributions_sc = zeros(Float32, length(residues_initial))
    sel_at_bb(at, r) = (at in r) & backbone(at) & selector(at)
    sel_at_sc(at, r) = (at in r) & sidechain(at) & selector(at)
    for iresidue in eachindex(residues_initial)
        rinit = residues_initial[iresidue]
        rfinal = residues_final[iresidue]
        rtype = threeletter(resname(rinit)) # convert non-standard residue names in types (e. g. HSD -> HIS)
        if !(haskey(protein_residues, rtype))
            throw(ArgumentError("""\n
                Found non-protein residue ($(resname(rinit))) in the selected atoms of SASA calculation.
                m-value calculations are only defined for protein residues.

            """))
        end
        if rtype != threeletter(resname(rfinal))
            throw(ArgumentError("""\n
                The residues of the initial and final state must be of the same type.
                Found residues at position $iresidue: $(resname(rinit)) and $(resname(rfinal))

            """))
        end
        bb_initial = sasa(sasa_initial, at -> sel_at_bb(at, rinit))
        sc_initial = sasa(sasa_initial, at -> sel_at_sc(at, rinit))
        bb_final = sasa(sasa_final, at -> sel_at_bb(at, rfinal))
        sc_final = sasa(sasa_final, at -> sel_at_sc(at, rfinal))
        bb_type, sc_type = tfe_asa(model, cosolvent, rtype)
        residue_contributions_bb[iresidue] = bb_type * (bb_final - bb_initial)
        residue_contributions_sc[iresidue] = sc_type * (sc_final - sc_initial)
    end
    bb = sum(residue_contributions_bb)
    sc = sum(residue_contributions_sc)
    tot = bb + sc
    return MValue{model}(length(residues_initial), tot, bb, sc, residue_contributions_bb, residue_contributions_sc)
end

#=
    tfe_asa(::Type{Urea}, restype::AbstractString)

Returns the transfer free energy per unit area (kcal/nm^2) for side chain and backbone
for a given amino acid type in urea solution according to the Tanford transfer model,
as implemente by Moeser and Horinek (https://pubs.acs.org/doi/10.1021/jp409934q#_i14).

=#
function tfe_asa(
    model::Type{<:MvalueModel},
    cosolvent::String,
    restype::AbstractString;
)
    col = cosolvent_column[cosolvent]
    if model == MoeserHorinek
        # united model: all bb ASA contributions are the same
        bb_contribution = tfe_sc_bb_moeser_and_horinek["BB"][col] / first(isolated_ASA["GLY"])
        sc_contribution = if restype == "GLY"
            0.0
        else
            tfe_sc_bb_moeser_and_horinek[restype][col] / last(isolated_ASA[restype])
        end
    elseif model == AutonBolen
        bb_contribution = tfe_sc_bb_auton_and_bolen["BB"][col] / first(isolated_ASA[restype])
        sc_contribution = if restype == "GLY"
            0.0
        else
            sc_contribution = tfe_sc_bb_auton_and_bolen[restype][col] / last(isolated_ASA[restype])
        end
    end
    # convert to kcal / nm^2 and return
    return bb_contribution / 1000, sc_contribution / 1000
end

include("./mvalue_exogenous.jl")
include("./testing.jl")