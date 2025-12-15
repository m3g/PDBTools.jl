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
    transfer_free_energy(atoms::AbstractVector{<:PDBTools.Atom}, cosolvent::String; kargs...)
    transfer_free_energy(sasa_atoms::SASA{3,<:AbstractVector{<:Atom}}, cosolvent::String; kargs...)

Calculates the transfer free energy (in 1M solution, in `kcal/mol`) using the Tanford transfer model,
as implemented by Moeser and Horinek [1] or by Auton and Bolen [2,3].

# Positional Arguments

- `atoms`:: Atoms of the system (a vector of PDBTools.Atom objects)

or

- `sasa_atoms::SASA{3,<:AbstractVector{<:Atom}}`: SASA object with solvent accessible surface area of the atoms.

and

- `cosolvent::String`: The cosolvent to consider. One of: $(join('"' .* sort!(unique(keys(PDBTools.cosolvent_column)) .* '"'; by=lowercase),", ")) (case insensitive).

# Keyword Arguments (optional)

- `model::Type{<:MValueModel}=AutonBolen`: The model to use for the calculation. Either `MoeserHorinek` or `AutonBolen`.
- `sel::Union{String,Function}=all`: Selection of atoms to consider in the calculation. Can be a selection string or a function that takes an `Atom` and returns a `Bool`.
- `backbone::Function = PDBTools.isbackbone`: Function to identify backbone atoms.
- `sidechain::Function = PDBTools.issidechain`: Function to identify side chain atoms.
- `parallel:Bool = true`: Set parallelization, requires starting Julia multithreaded.

# Returns

A `TransferFreeEnergy` object, with fields:

- `ntatoms::Int`: Number of atoms considered.
- `tot::Float32`: Total m-value (kcal/mol/M).
- `bb::Float32`: Backbone contribution to the m-value (kcal/mol/M).
- `sc::Float32`: Side chain contribution to the m-value (kcal/mol/M).
- `residue_contributions_bb::Vector{Float32}`: Backbone contributions of each residue to the m-value.
- `residue_contributions_sc::Vector{Float32}`: Side-chain contributions of each residue to the m-value.
- `cosolvent::STring`: The cosolvent considered.

## Example

```julia
using PDBTools
# directly from the atoms:
prot = read_pdb("native.pdb")
transfer_free_energy(prot, "urea")
# or providing the SASA object:
sasa_prot = sasa_particles(prot)
transfer_free_energy(sasa_prot, "urea")
```

!!! warning
    When providing the vector of atoms, the SASA calculation will be performed using 
    Creamer united atom types. On the other hand, if the SASAs are computed previously,
    they will correspond to the method used which, by default, consider all atoms of 
    the structure, including hydrogens, and uses standard vdW atoms.

## References

1. https://doi.org/10.1021/jp409934q
2. https://doi.org/10.1016/s0076-6879(07)28023-1
3. https://www.pnas.org/doi/10.1073/pnas.0706251104

"""
function transfer_free_energy(
    atoms::AbstractVector{<:Atom},
    cosolvent::String;
    parallel=true,
    kargs...
)
    sasa_ats = sasa_particles(
        atoms;
        atom_type = creamer_atom_type,
        atom_radius_from_type = at -> creamer_atomic_radii[at],
        parallel
    ) 
    return transfer_free_energy(sasa_ats, cosolvent; parallel, kargs...)
end

function transfer_free_energy(
    sasa_ats::SASA{3,<:AbstractVector{<:Atom}},
    cosolvent::String;
    sel::Union{String,Function}=all,
    model::Type{<:MValueModel}=AutonBolen,
    backbone::F1=isbackbone,
    sidechain::F2=issidechain,
    parallel::Bool=true,
) where {F1<:Function,F2<:Function}
    selector = Select(sel)
    ats = select(sasa_ats.particles, selector)
    residues = collect(eachresidue(ats))
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
            if !(haskey(protein_residues, rtype))
                throw(ArgumentError("""\n
                    Found non-protein residue ($(resname(res))) in the selected atoms of SASA calculation.
                    m-value calculations are only defined for protein residues.
    
                """))
            end
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