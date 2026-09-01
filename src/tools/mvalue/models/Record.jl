#=

Record model

https://www.pnas.org/doi/abs/10.1073/pnas.1109372108?doi=10.1073%2Fpnas.1109372108

=#

export MTRecord
struct MTRecord <: MValueModel end

export MTRecordDenaturedModel

# Do not user underscores (_) in the following names:
const cosolvent_column_MTRecord = OrderedDict(
    "betaine" => 1,
    "urea" => 2,
    "tmao" => 3,
    "proline" => 4,
    "trehalose" => 5,
    "tetraeg" => 6,
    "glycerol" => 7,
)

#=

Each ASA is divided into contributions from
seven coarse-grained surface types (Tables S2 and S3): aliphatic
carbon, aromatic carbon, hydroxyl oxygen, amide oxygen,
carboxylate oxygen, amide nitrogen, and cationic nitrogen. All
amino acid surfaces appear well-described by these surface types
except His ring N, which is assigned as cationic N, Trp ring N,
which is assigned as amide N, and Met and Cys S, which are ig-
nored because the ΔASA of sulfur in unfolding the proteins in
Table S3 is negligibly small (0–2% of the total ΔASA).

=#
const _record_aromatic_carbons = Set{Tuple{String,String}}([
    ("PHE", "CG"), ("PHE", "CD1"), ("PHE", "CD2"), ("PHE", "CE1"), ("PHE", "CE2"), ("PHE", "CZ"),
    ("TYR", "CG"), ("TYR", "CD1"), ("TYR", "CD2"), ("TYR", "CE1"), ("TYR", "CE2"), ("TYR", "CZ"),
    ("TRP", "CG"), ("TRP", "CD1"), ("TRP", "CD2"), ("TRP", "CE2"), ("TRP", "CE3"), ("TRP", "CZ2"), ("TRP", "CZ3"), ("TRP", "CH2"),
    ("HIS", "CG"), ("HIS", "CE1"), ("HIS", "CD2"),
])

const _record_hydroxyl_oxygens = Set{Tuple{String,String}}([
    ("SER", "OG"),
    ("THR", "OG1"),
    ("TYR", "OH"),
])

const _record_amide_oxygens = Set{Tuple{String,String}}([
    ("ASN", "OD1"),
    ("GLN", "OE1"),
])

const _record_carboxylate_oxygens = Set{Tuple{String,String}}([
    ("ASP", "OD1"),
    ("ASP", "OD2"),
    ("GLU", "OE1"),
    ("GLU", "OE2"),
])

const _record_amide_nitrogens = Set{Tuple{String,String}}([
    ("ASN", "ND2"),
    ("GLN", "NE2"),
    ("TRP", "NE1"), # Explicit exception in the Record model.
])

const _record_cationic_nitrogens = Set{Tuple{String,String}}([
    ("LYS", "NZ"),
    ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
    ("HIS", "ND1"), ("HIS", "NE2"), # Explicit exception in the Record model.
])

#=
    record_surface_type(atom::Atom; nterminal::Bool=false, cterminal::Bool=false)

Map a protein atom to one of the seven coarse-grained surface types from the
Record transfer model:

- `:aliphatic_carbon`
- `:aromatic_carbon`
- `:hydroxyl_oxygen`
- `:amide_oxygen`
- `:carboxylate_oxygen`
- `:amide_nitrogen`
- `:cationic_nitrogen`

Returns `nothing` for sulfur and hydrogen atoms, which are ignored by this model.

Since this classification is otherwise a pure function of the atom itself, terminal
groups need residue-position context passed in explicitly: set `nterminal=true` for
the backbone "N" of a chain's first residue (free α-amino group, protonated and
cationic rather than an amide nitrogen), and `cterminal=true` for the backbone "O" of
a residue whose terminal carboxylate is present (i.e. it also carries an OXT/OT1/OT2
atom), so that "O" is classified alongside OXT as `:carboxylate_oxygen` instead of
`:amide_oxygen`. Both default to `false`, which reproduces the original, purely
per-atom, classification (used e.g. by the `record_*` macro keywords, which have no
residue context available).
=#
function record_surface_type(at::Atom; nterminal::Bool=false, cterminal::Bool=false)
    restype = String(uppercase(threeletter(StringType, resname(at))))
    atname = String(uppercase(name(at)))
    elem = element(at)
    isnothing(elem) && throw(ArgumentError("Could not determine element of atom $(name(at)) in residue $(resname(at))."))
    el = String(uppercase(elem))
    key = (restype, atname)

    if el == "C"
        return key in _record_aromatic_carbons ? :aromatic_carbon : :aliphatic_carbon
    end

    if el == "O"
        if atname in ("OXT", "OT1", "OT2") || key in _record_carboxylate_oxygens
            return :carboxylate_oxygen
        elseif key in _record_hydroxyl_oxygens
            return :hydroxyl_oxygen
        elseif atname == "O"
            return cterminal ? :carboxylate_oxygen : :amide_oxygen
        elseif key in _record_amide_oxygens
            return :amide_oxygen
        end
        return :amide_oxygen
    end

    if el == "N"
        if key in _record_cationic_nitrogens
            return :cationic_nitrogen
        elseif atname == "N"
            return nterminal ? :cationic_nitrogen : :amide_nitrogen
        elseif key in _record_amide_nitrogens
            return :amide_nitrogen
        end
        return :amide_nitrogen
    end

    if el == "S" || el == "H"
        return nothing
    end

    throw(ArgumentError("Atom $(name(at)) of residue $(resname(at)) has unsupported element $el for Record model."))
end

const _record_macro_keyword_functions = (
    "record_aliphatic_carbon" => (at -> record_surface_type(at) === :aliphatic_carbon),
    "record_aromatic_carbon" => (at -> record_surface_type(at) === :aromatic_carbon),
    "record_hydroxyl_oxygen" => (at -> record_surface_type(at) === :hydroxyl_oxygen),
    "record_amide_oxygen" => (at -> record_surface_type(at) === :amide_oxygen),
    "record_carboxylate_oxygen" => (at -> record_surface_type(at) === :carboxylate_oxygen),
    "record_amide_nitrogen" => (at -> record_surface_type(at) === :amide_nitrogen),
    "record_cationic_nitrogen" => (at -> record_surface_type(at) === :cationic_nitrogen),
)

function _register_record_macro_keywords!()
    for (kname, getter) in _record_macro_keyword_functions
        if !any(k -> k.name == kname, macro_keywords)
            push!(macro_keywords, MacroKeyword(kname, getter))
        end
    end
    return nothing
end

_register_record_macro_keywords!()

const _record_alpha_surface = Dict{String,Dict{Symbol,Float32}}(
    # Values from Guinn et al. (PNAS 2011), Table 1.
    # Stored as 10^4 * alpha_i (units m^-1 A^-2).
    "urea" => Dict{Symbol,Float32}(
        :aromatic_carbon => -8.9f0,
        :amide_oxygen => -8.7f0,
        :carboxylate_oxygen => -4.0f0,
        :amide_nitrogen => -3.2f0,
        :hydroxyl_oxygen => -2.5f0,
        :aliphatic_carbon => -1.1f0,
        :cationic_nitrogen => 1.8f0,
    ),
    "betaine" => Dict{Symbol,Float32}(
        :aromatic_carbon => -23.0f0,
        :amide_oxygen => 28.0f0,
        :carboxylate_oxygen => 29.0f0,
        :amide_nitrogen => -20.0f0,
        :hydroxyl_oxygen => 1.0f0,
        :aliphatic_carbon => 3.0f0,
        :cationic_nitrogen => -12.0f0,
    ),
    # TMAO data from: https://doi.org/10.1016/j.bpj.2016.09.035
    "tmao" => Dict{Symbol,Float32}( 
        :aromatic_carbon => 22.6f0,
        :amide_oxygen => 9.3f0,
        :carboxylate_oxygen => 52.3f0,
        :amide_nitrogen => 11.6f0,
        :hydroxyl_oxygen => 9.9f0,
        :aliphatic_carbon => 17.8f0,
        :cationic_nitrogen => -4.6f0,
    ),
    # Proline data from: https://doi.org/10.1021/bi400683y
    "proline" => Dict{Symbol,Float32}( 
        :aromatic_carbon => -9.2f0,
        :amide_oxygen => 14.5f0,
        :carboxylate_oxygen => 16.6f0,
        :amide_nitrogen => -11.8f0,
        :hydroxyl_oxygen => -0.7f0,
        :aliphatic_carbon => 5.3f0,
        :cationic_nitrogen => -12.6f0,
    ),
    # Trehalose data from: https://doi.org/10.1016/j.bpj.2015.05.037
    "trehalose" => Dict{Symbol,Float32}(
        :aromatic_carbon => 5.9,
        :amide_oxygen => -19.6,
        :carboxylate_oxygen => -28.2,
        :amide_nitrogen => -4.7,
        :hydroxyl_oxygen => -0.8,
        :aliphatic_carbon => 22.4,
        :cationic_nitrogen => 12.9,
    ),
)

#
# TetraEG (tetraethylene glycol) and glycerol group interaction potentials, from
# Table 1 ("Group Interaction Potentials (α_solute,i, cal mol⁻¹ molal⁻¹ Å⁻²)") of the
# paper reporting interactions of TetraEG, glycerol, and end/interior groups of PEG with
# protein groups and salt ions. Unlike `_record_alpha_surface` above (Guinn et al.
# convention, stored as 10^4 * alpha_i and rescaled by R*T in `model_combination_rule`),
# these values are already given directly in cal mol⁻¹ molal⁻¹ Å⁻², so they are used as-is.
# The table's "carboxylic acid O" row has no counterpart in `record_surface_type` (which,
# following Guinn et al., only has a single, deprotonated `:carboxylate_oxygen` class) and
# is therefore omitted here.
#
# Data from: https://doi.org/10.1021/acs.biochem.5b00246
#
const _record_alpha_surface_cal = Dict{String,Dict{Symbol,Float32}}(
    "tetraeg" => Dict{Symbol,Float32}(
        :aliphatic_carbon => -0.349f0,
        :aromatic_carbon => -2.66f0,
        :hydroxyl_oxygen => 0.843f0,
        :amide_oxygen => 2.94f0,
        :carboxylate_oxygen => 3.59f0,
        :amide_nitrogen => -1.67f0,
        :cationic_nitrogen => -0.805f0,
    ),
    "glycerol" => Dict{Symbol,Float32}(
        :aliphatic_carbon => 0.0548f0,
        :aromatic_carbon => -0.431f0,
        :hydroxyl_oxygen => 0.0305f0,
        :amide_oxygen => 0.826f0,
        :carboxylate_oxygen => 0.467f0,
        :amide_nitrogen => -0.491f0,
        :cationic_nitrogen => -0.245f0,
    ),
)

const _record_R_cal = 1.9872041f0
const _record_default_T = 298.15f0

#=
    model_combination_rule(::Type{MTRecord}, cosolvent, surface_type::Symbol)

Return the Record surface interaction potential contribution in `cal mol⁻¹ A⁻²` for
a coarse-grained surface type, computed from Guinn et al. Table 1 and Eq. 4 (for
`_record_alpha_surface` cosolvents), or read directly from `_record_alpha_surface_cal`
(TetraEG and glycerol, already in cal mol⁻¹ molal⁻¹ Å⁻²).
=#
function model_combination_rule(::Type{MTRecord}, cosolvent, surface_type::Symbol)
    cos = lowercase(cosolvent)
    if haskey(_record_alpha_surface_cal, cos)
        α = get(_record_alpha_surface_cal[cos], surface_type, nothing)
        isnothing(α) && throw(ArgumentError("Unsupported MTRecord surface type: $surface_type"))
        return α
    end
    haskey(_record_alpha_surface, cos) || throw(ArgumentError("Unsupported cosolvent for MTRecord: $cosolvent"))
    α_1e4 = get(_record_alpha_surface[cos], surface_type, nothing)
    isnothing(α_1e4) && throw(ArgumentError("Unsupported MTRecord surface type: $surface_type"))
    α = α_1e4 * 1f-4
    return _record_R_cal * _record_default_T * α
end

"""
    MTRecordDenaturedModel

Type that specifies that a m-value calculation will consider the `MTRecord` transfer
model (Guinn et al., PNAS 2011). Unlike the other transfer models, `MTRecord` was
parameterized against, and is meant to be used with, an explicit atomistic model of the
unfolded state: this type wraps a protein's native (folded) structure together with its
fully-extended (all-trans, phi = psi = 180°, see `extended_chain`) chain, and the
atom-resolved SASAs of both, precomputed once at construction. This type is used as the
first input variable of the `mvalue` function.

Construction:

```
MTRecordDenaturedModel(atoms::AbstractVector{<:Atom})
```

builds the extended chain internally (via `extended_chain`) and computes both SASAs.

Use the `MTRecordDenaturedModel` model as the first input argument of `mvalue`, for example:

```
mvalue(MTRecordDenaturedModel(prot), "urea")
```

to obtain the estimated m-value of denaturation in "urea".

Reference:

Guinn EJ, Pegram LM, Capp MW, Pollock MN, Record MT Jr. **Quantifying why urea is a protein
denaturant, whereas glycine betaine is a protein stabilizer.** *PNAS.* 2011;108:16932-16937.
doi: 10.1073/pnas.1109372108.

"""
struct MTRecordDenaturedModel{T1,T2,S1,S2}
    native_chain::T1
    extended_chain::T2
    sasa_native::S1
    sasa_ext::S2
end

function MTRecordDenaturedModel(p::AbstractVector{<:Atom})
    p_ext = extended_chain(p)
    sasa_native = sasa_particles(RichardsUnitedAtomRadii, p)
    sasa_ext = sasa_particles(RichardsUnitedAtomRadii, p_ext)
    return MTRecordDenaturedModel(p, p_ext, sasa_native, sasa_ext)
end

function Base.show(io::IO, m::MTRecordDenaturedModel)
    print(io, "MTRecordDenaturedModel wrapping a $(length(m.native_chain))-atom native/extended chain pair")
end

"""
    mvalue(m::MTRecordDenaturedModel, cosolvent::AbstractString; alpha=1.0)

Compute Record-model m-values from coarse-grained surface interaction potentials, using
Eq. 4 (Guinn et al. 2011) and the atom-resolved ΔASA between the native and fully-extended
chains wrapped in `m` (that is, ASA of the extended chain minus ASA of the native chain, for
each atom), following the original Record model.

The optional `alpha` keyword scales the extended-chain (denatured-state) ASA of each atom
before taking the difference with the native-state ASA. The default, `alpha=1.0`, uses the
extended-chain ASA as computed, with no adjustment; this reproduces the urea benchmarks of
Guinn et al. (Table S3) well. There is no equivalent extended-chain validation set for
betaine in the paper, so `alpha` is exposed to allow empirically fitting the denatured-state
ASA scale to experimental betaine data, if needed.

"""
function mvalue(m::MTRecordDenaturedModel, cosolvent::AbstractString; alpha=1.0, kargs...)
    tfe_n = transfer_free_energy(MTRecord, m.sasa_native, cosolvent; kargs...)
    tfe_d = transfer_free_energy(MTRecord, m.sasa_ext, cosolvent; kargs...)
    mval = alpha * tfe_d - tfe_n
    return MValue{MTRecord}(
        mval.nresidues, mval.tot, mval.bb, mval.sc,
        mval.residue_contributions_bb, mval.residue_contributions_sc,
        mval.cosolvent
    )
end

"""
    transfer_free_energy(::Type{MTRecord}, atoms::AbstractVector{<:Atom}, cosolvent::AbstractString; kargs...)

Compute transfer free energies from a fixed (native) structure using the Record
surface-type model and atom-resolved SASAs.

"""
function transfer_free_energy(
    ::Type{MTRecord},
    atoms::AbstractVector{<:Atom},
    cosolvent::AbstractString;
    backbone::F1=isbackbone,
    sel::Union{String,Function}=all,
    sidechain::F2=issidechain,
    parallel::Bool=true,
    unitcell=nothing,
) where {F1<:Function,F2<:Function}
    sasa_ats = sasa_particles(RichardsUnitedAtomRadii, atoms; unitcell)
    return transfer_free_energy(
        MTRecord,
        sasa_ats,
        cosolvent;
        backbone,
        sel,
        sidechain,
        parallel,
    )
end

"""
    transfer_free_energy(::Type{MTRecord}, sasa_ats::SASA{RichardsUnitedAtomRadii}, cosolvent::AbstractString; kargs...)

Compute fixed-state transfer free energies using Record interaction potentials
and atom-resolved coarse-grained surface classes.

"""
function transfer_free_energy(
    ::Type{MTRecord},
    sasa_ats::SASA{RichardsUnitedAtomRadii},
    cosolvent::AbstractString;
    backbone::F1=isbackbone,
    sel::Union{String,Function}=all,
    sidechain::F2=issidechain,
    parallel::Bool=true,
) where {F1<:Function,F2<:Function}
    selector = Select(sel)
    residues = collect(eachresidue(select(sasa_ats.particles, selector)))
    cos = lowercase(cosolvent)
    haskey(cosolvent_column_MTRecord, cos) || throw(ArgumentError("Unsupported cosolvent for MTRecord: $cosolvent"))

    residue_contributions_bb = zeros(Float32, length(residues))
    residue_contributions_sc = zeros(Float32, length(residues))
    nchunks = parallel ? Threads.nthreads() : 1
    @sync for inds_chunk in ChunkSplitters.index_chunks(residues; n=nchunks)
        Threads.@spawn for iresidue in inds_chunk
            res = residues[iresidue]
            bb_atoms = [at for at in res if backbone(at) && selector(at)]
            sc_atoms = [at for at in res if sidechain(at) && selector(at)]
            bb_contrib = 0.0f0
            sc_contrib = 0.0f0

            # Terminal groups need residue-position context: the first residue of each
            # chain carries a free (cationic) N-terminal amine, and a residue carrying an
            # OXT/OT1/OT2 atom has a C-terminal carboxylate, so its "O" is chemically
            # equivalent to OXT rather than an amide oxygen.
            is_nterminal = iresidue == 1 || chain(res) != chain(residues[iresidue-1]) || model(res) != model(residues[iresidue-1])
            is_cterminal = any(at -> String(uppercase(name(at))) in ("OXT", "OT1", "OT2"), res)

            for at in bb_atoms
                stype = record_surface_type(at; nterminal=is_nterminal, cterminal=is_cterminal)
                isnothing(stype) && continue
                σ = model_combination_rule(MTRecord, cos, stype)
                bb_contrib += σ * sasa(sasa_ats, x -> x === at)
            end

            for at in sc_atoms
                stype = record_surface_type(at)
                isnothing(stype) && continue
                σ = model_combination_rule(MTRecord, cos, stype)
                sc_contrib += σ * sasa(sasa_ats, x -> x === at)
            end

            residue_contributions_bb[iresidue] = bb_contrib / 1000f0
            residue_contributions_sc[iresidue] = sc_contrib / 1000f0
        end
    end

    bb = sum(residue_contributions_bb)
    sc = sum(residue_contributions_sc)
    tot = bb + sc
    return TransferFreeEnergy{MTRecord}(length(residues), tot, bb, sc, residue_contributions_bb, residue_contributions_sc, cos)
end

"""
    transfer_free_energy(sasa_ats::SASA{RichardsUnitedAtomRadii}, cosolvent::AbstractString; model=MTRecord, kargs...)

Compute transfer free energies from a precomputed SASA using Richards united atom radii.
This is the `RichardsUnitedAtomRadii` counterpart of
`transfer_free_energy(sasa_ats::SASA{CreamerUnitedAtomRadii}, ...)`: since Richards radii are
only meaningful for the `MTRecord` model, `model` must be `MTRecord` (its default).

"""
function transfer_free_energy(
    sasa_ats::SASA{RichardsUnitedAtomRadii},
    cosolvent::AbstractString;
    model::Type{<:MValueModel}=MTRecord,
    kargs...
)
    if model != MTRecord
        throw(ArgumentError("""\n
            SASA computed with RichardsUnitedAtomRadii can only be used with the MTRecord model.
            Got model: $(modelname(model))

        """))
    end
    return transfer_free_energy(MTRecord, sasa_ats, cosolvent; kargs...)
end

