export mvalue
export MoeserHorinek, AutonBolen

abstract type MvalueModel end
struct MoeserHorinek <: MvalueModel end
struct AutonBolen <: MvalueModel end

include("./data.jl")

"""
    mvalue(; model=MoeserHorinek, cosolvent="urea", atoms:AbstractVector{<:PDBTools.Atom}, sasas, type=1)

Calculates the m-value (transfer free energy of a protein in 1M solution) using the Tanford transfer model,
as implemented by Moeser and Horinek [1] or by Auton and Bolen [2,3].

# Arguments

- `model`: The model to be used. Must be `MoeserHorinek` or `AutonBolen`. `MoeserHorinek` is only implemented for `cosolvent="urea"`,
   and should be more precise in that case. Other solvents are available for `AutonBolen`.
- `cosolvent::String`: One of $(join('"' .* sort!(unique(keys(PDBTools.cosolvent_column)) .* '"'; by=lowercase),", "))
- `atoms::AbstractVector{<:PDBTools.Atom}`: Vector containing the atoms of the structure.
- `sasas::Dict{String, Dict{Symbol, Float64}}`: A dictionary containing the change in solvent accessible surface area (SASA)
  upon denaturation for each amino acid type. This data can be obtained from the m-value server or calculated using GROMACS:
    - The output of the server can be parsed using the `parse_mvalue_server_sasa` function defined in this module.
    - Compute the SASA with `delta_sasa_per_restype`, a SASA calculation utility implemented in PDBTools.jl.
    - SASA values can be calculated using GROMACS with the `gmx_delta_sasa_per_restype` function defined in this module.

- `type::Int`: Specifies which SASA value to use from the provided data, because the server provides minimum, average, and maximum values,
    according to different denatured models for the protein. The recommended value is `2` for comparison with experimental data.
    Normally, GROMACS calculations will provide a single value, so `type=1` should be used in that case.

# Returns

A named tuple with the following fields:
- `tot`: Total transfer free energy (kcal/mol).
- `bb`: Contribution from the backbone (kcal/mol).
- `sc`: Contribution from the side chains (kcal/mol).
- `restype`: A dictionary with the transfer free energy contributions per residue type (kcal/mol).

Each entry in the dictionary is a named tuple with `bb` and `sc` fields representing the backbone and side chain contributions, respectively.

# Example calls

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
    selection::Union{String,Function}=all;
    model::Type{<:MvalueModel}=MoeserHorinek,
    cosolvent::String="urea",
)
    if selection isa String
        sel_func = Select(selection)
    else
        sel_func = selection
    end
    ats_initial = select(sasa_initial.particles, sel_func)
    ats_final = select(sasa_final.particles, sel_func)
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
        bb_initial = sasa(sasa_initial, at -> (resname(at) == rtype) & isbackbone(at) & sel_func(at) )
        sc_initial = sasa(sasa_initial, at -> (resname(at) == rtype) & issidechain(at) & sel_func(at) )
        bb_final = sasa(sasa_final, at -> (resname(at) == rtype) & isbackbone(at) & sel_func(at) )
        sc_final = sasa(sasa_final, at -> (resname(at) == rtype) & issidechain(at) & sel_func(at) )
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

