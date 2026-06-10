export transfer_free_energy

struct TransferFreeEnergy{T}
    nresidues::Int
    tot::Float32
    bb::Float32
    sc::Float32
    residue_contributions_bb::Vector{Float32}
    residue_contributions_sc::Vector{Float32}
    cosolvent::String
end
function Base.show(io::IO, ::MIME"text/plain", t::TransferFreeEnergy)
    print(io, chomp("""
    $(typeof(t)) - $(t.nresidues) residues to 1M \"$(t.cosolvent)\".
        Total transfer free energy: $(t.tot) kcal mol⁻¹
        Backbone contributions: $(t.bb) kcal mol⁻¹
        Side-chain contributions: $(t.sc) kcal mol⁻¹
    """))
end

"""
    transfer_free_energy(atoms::AbstractVector{<:PDBTools.Atom}, cosolvent::AbstractString; kargs...)

Calculates the transfer free energy (in 1M solution, in `kcal/mol`) using the Tanford transfer model,
as implemented by Moeser and Horinek [1] or by Auton and Bolen [2,3].

# Positional Arguments

- `atoms`:: Atoms of the system (a vector of PDBTools.Atom objects)
- `cosolvent::AbstractString`: The cosolvent to consider. 

$(_available_cosolvents())

# Keyword Arguments (optional)

- `model::Type{<:MValueModel}=AutonBolen`: The model to use for the calculation. 
- `sel::Union{String,Function}=all`: Selection of atoms to consider in the calculation. Can be a selection string or a function that takes an `Atom` and returns a `Bool`.
- `backbone::Function = PDBTools.isbackbone`: Function to identify backbone atoms.
- `sidechain::Function = PDBTools.issidechain`: Function to identify side chain atoms.
- `parallel:Bool = true`: Set parallelization, requires starting Julia multithreaded.
- `unitcell=nothing`: if periodic boundary conditions are used, provide a 3x3 matrix with
  the unitcell, or alternatively a vector of length 3 with the sides, for orthorhombic cells.

# Returns

A `TransferFreeEnergy` object, with fields:

- `ntatoms::Int`: Number of atoms considered.
- `tot::Float32`: Total m-value (kcal/mol/M).
- `bb::Float32`: Backbone contribution to the m-value (kcal/mol/M).
- `sc::Float32`: Side chain contribution to the m-value (kcal/mol/M).
- `residue_contributions_bb::Vector{Float32}`: Backbone contributions of each residue to the m-value.
- `residue_contributions_sc::Vector{Float32}`: Side-chain contributions of each residue to the m-value.
- `cosolvent::String`: The cosolvent considered.

## Example

```julia
using PDBTools
prot = read_pdb("native.pdb")
transfer_free_energy(prot, "urea")
```
## References

1. https://doi.org/10.1021/jp409934q
2. https://doi.org/10.1016/s0076-6879(07)28023-1
3. https://www.pnas.org/doi/10.1073/pnas.0706251104

"""
function transfer_free_energy(
    atoms::AbstractVector{<:Atom},
    cosolvent::AbstractString;
    model::Type{<:MValueModel}=AutonBolen,
    backbone::F1=isbackbone,
    sel::Union{String,Function}=all,
    sidechain::F2=issidechain,
    parallel::Bool=true,
    unitcell=nothing,
) where {F1<:Function,F2<:Function}
    sasa_ats = sasa_particles(CreamerUnitedAtomRadii, atoms; unitcell)
    return transfer_free_energy(
        sasa_ats, cosolvent;
        model, backbone, sel, sidechain, parallel
    )
end

function transfer_free_energy(sasa_ats::SASA{T1}, cosolvent; kargs...) where {T1}
    throw(ArgumentError("""\n
        To computie m-values or transfer free energies the SASA computation 
        must use CreamerUnitedAtomRadii. For example, use:

            s = sasa_particles(CreamerUnitedAtomRadii, atoms)

        Got atomic radii type: $T1

    """))
end

"""
    transfer_free_energy(
        sasa_ats::SASA{CreamerUnitedAtomRadii},
        cosolvent::AbstractString;
        model::Type{<:MValueModel}=AutonBolen,
        backbone::F1=isbackbone,
        sel::Union{String,Function}=all,
        sidechain::F2=issidechain,
        parallel::Bool=true,
    )
    

Compute transfer free energies from precomputed solvent accessible surface areas. The SASAs must
have been computed with `CreamerUnitedAtomRadii` radii.

"""
function transfer_free_energy(
    sasa_ats::SASA{CreamerUnitedAtomRadii},
    cosolvent::AbstractString;
    model::Type{<:MValueModel}=AutonBolen,
    backbone::F1=isbackbone,
    sel::Union{String,Function}=all,
    sidechain::F2=issidechain,
    parallel::Bool=true,
) where {F1,F2}
    selector = Select(sel)
    residues = collect(eachresidue(select(sasa_ats.particles, selector)))
    cosolvent = lowercase(cosolvent)
    residue_contributions_bb = zeros(Float32, length(residues))
    residue_contributions_sc = zeros(Float32, length(residues))
    nchunks = parallel ? Threads.nthreads() : 1
    @sync for inds_chunk in ChunkSplitters.index_chunks(residues; n=nchunks)
        Threads.@spawn for iresidue in inds_chunk
            sel_at_bb = _AtomSelector(at -> backbone(at) & selector(at))
            sel_at_sc = _AtomSelector(at -> sidechain(at) & selector(at))
            res = residues[iresidue]
            rtype = threeletter(resname(res)) # convert non-standard residue names in types (e. g. HSD -> HIS)
            sel_at_sc.residue[] = res
            sel_at_bb.residue[] = res
            bb_res = sasa(sasa_ats, sel_at_bb) 
            sc_res = sasa(sasa_ats, sel_at_sc)
            bb_type, sc_type = tfe_asa(model, cosolvent, rtype)
            residue_contributions_bb[iresidue] = bb_type * bb_res
            residue_contributions_sc[iresidue] = sc_type * sc_res
        end
    end
    bb = sum(residue_contributions_bb)
    sc = sum(residue_contributions_sc)
    tot = bb + sc
    return TransferFreeEnergy{model}(length(residues), tot, bb, sc, residue_contributions_bb, residue_contributions_sc, cosolvent)
end