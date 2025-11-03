export mvalue
export MoeserHorinek, AutonBolen

abstract type MvalueModel end
struct MoeserHorinek <: MvalueModel end
struct AutonBolen <: MvalueModel end

include("./data.jl")

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

Calculates the m-value (transfer free energy of a protein in 1M solution) using the Tanford transfer model,
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

A named tuple with the following fields:

- `tot::Float64`: Total m-value (kcal/mol/M).
- `bb::Float64`: Backbone contribution to the m-value (kcal/mol/M).
- `sc::Float64`: Side chain contribution to the m-value (kcal/mol/M).
- `restype::Dict{String,@NamedTuple{bb::Float64, sc::Float64}}`: Dictionary with per-residue type contributions to the m-value.

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
    backbone::Function = isbackbone,
    sidechain::Function = issidechain,
)
    selector = Select(sel)
    ats_initial = filter(selector, sasa_initial.particles)
    ats_final = filter(selector, sasa_final.particles)
    residue_names = unique(vcat(resname.(ats_initial), resname.(ats_final)))

    DeltaG_per_residue = Dict{String,@NamedTuple{bb::Float64, sc::Float64}}()
    for rname in residue_names
        rtype = threeletter(rname) # convert non-standard residue names in types (e. g. HSD -> HIS)
        if !(haskey(protein_residues, rtype))
            throw(ArgumentError("""\n
                Found non-protein residue ($rname) in the selected atoms of SASA calculation.
                m-value calculations are only defined for protein residues.
            
            """))
        end
        bb_initial = sasa(sasa_initial, at -> (resname(at) == rtype) & backbone(at) & selector(at) )
        sc_initial = sasa(sasa_initial, at -> (resname(at) == rtype) & sidechain(at) & selector(at) )
        bb_final = sasa(sasa_final, at -> (resname(at) == rtype) & backbone(at) & selector(at) )
        sc_final = sasa(sasa_final, at -> (resname(at) == rtype) & sidechain(at) & selector(at) )
        DeltaG_per_residue[rtype] = (
            bb=(last(tfe_asa(model, cosolvent, rtype)) * (bb_final - bb_initial) / 100),
            sc=(first(tfe_asa(model, cosolvent, rtype)) * (sc_final - sc_initial) / 100),
        )
    end
    DeltaG_BB = sum(getfield(DeltaG_per_residue[key], :bb) for key in keys(DeltaG_per_residue); init=0.f0)
    DeltaG_SC = sum(getfield(DeltaG_per_residue[key], :sc) for key in keys(DeltaG_per_residue); init=0.f0)
    DeltaG = DeltaG_BB + DeltaG_SC
    return (tot=DeltaG, bb=DeltaG_BB, sc=DeltaG_SC, restype=DeltaG_per_residue)
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
    else
        error("model must be either MoeserHorinek or AutonBolen")
    end
    # convert to kcal / nm^2 and return
    return sc_contribution / 10, bb_contribution / 10
end

include("./mvalue_exogenous.jl")
include("./testing.jl")