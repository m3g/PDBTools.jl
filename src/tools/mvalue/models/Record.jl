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

"""
    record_surface_type(atom::Atom)

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
"""
function record_surface_type(at::Atom)
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
        elseif atname == "O" || key in _record_amide_oxygens
            return :amide_oxygen
        end
        return :amide_oxygen
    end

    if el == "N"
        if key in _record_cationic_nitrogens
            return :cationic_nitrogen
        elseif atname == "N" || key in _record_amide_nitrogens
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
)

const _record_R_cal = 1.9872041f0
const _record_default_T = 298.15f0

"""
    model_combination_rule(::Type{MTRecord}, cosolvent, surface_type::Symbol)

Return the Record surface interaction potential contribution in `cal mol⁻¹ A⁻²` for
a coarse-grained surface type, computed from Guinn et al. Table 1 and Eq. 4.
"""
function model_combination_rule(::Type{MTRecord}, cosolvent, surface_type::Symbol)
    cos = lowercase(cosolvent)
    haskey(_record_alpha_surface, cos) || throw(ArgumentError("Unsupported cosolvent for MTRecord: $cosolvent"))
    α_1e4 = get(_record_alpha_surface[cos], surface_type, nothing)
    isnothing(α_1e4) && throw(ArgumentError("Unsupported MTRecord surface type: $surface_type"))
    α = α_1e4 * 1f-4
    return _record_R_cal * _record_default_T * α
end

function model_combination_rule(::Type{MTRecord}, cosolvent, restype::AbstractString)
    throw(ArgumentError("""
        MTRecord is an atom/surface-type model and cannot be combined per-residue type.
        Use mvalue(MTRecordDenaturedModel(atoms[, type]), cosolvent) instead.
    """))
end

struct MTRecordDenaturedModel{C}
    creamer::C
end

function MTRecordDenaturedModel(c::CreamerDenaturedModel)
    return MTRecordDenaturedModel{typeof(c)}(c)
end

function MTRecordDenaturedModel(
    atoms::AbstractVector{<:Atom},
    type::Int=3;
    sasa_parameterization=:original,
)
    return MTRecordDenaturedModel(
        CreamerDenaturedModel(atoms, type; sasa_parameterization)
    )
end

function Base.show(io::IO, m::MTRecordDenaturedModel)
    print(io, "MTRecordDenaturedModel wrapping $(m.creamer)")
end

_record_denatured_sc_bb(c::CreamerDenaturedModel, rtype::String) = begin
    cr = _sasa_parameterization(c.sasa_parameterization)[rtype]
    bb = c.type == 1 ? cr.bb_lower : c.type == 2 ? 0.5f0 * (cr.bb_lower + cr.bb_upper) : cr.bb_upper
    sc = c.type == 1 ? cr.sc_lower : c.type == 2 ? 0.5f0 * (cr.sc_lower + cr.sc_upper) : cr.sc_upper
    return bb, sc
end

"""
    mvalue(m::MTRecordDenaturedModel, cosolvent::AbstractString)

Compute Record-model m-values from coarse-grained surface interaction potentials
using Eq. 4 (Guinn et al. 2011) and folded SASAs from the wrapped
`CreamerDenaturedModel`.
"""
function mvalue(m::MTRecordDenaturedModel, cosolvent::AbstractString; alpha=1.15)
    cos = lowercase(cosolvent)
    haskey(cosolvent_column_MTRecord, cos) || throw(ArgumentError("Unsupported cosolvent for MTRecord: $cosolvent"))

    c = m.creamer
    residues = collect(eachresidue(c.sasa_atoms.particles))
    n = length(residues)
    residue_contributions_bb = zeros(Float32, n)
    residue_contributions_sc = zeros(Float32, n)

    for (i, res) in enumerate(residues)
        rtype = threeletter(resname(res))
        bb_denatured, sc_denatured = alpha .* _record_denatured_sc_bb(c, rtype)

        bb_atoms = [at for at in res if isbackbone(at)]
        sc_atoms = [at for at in res if issidechain(at)]

        bb_folded = 0.0f0
        sc_folded = 0.0f0
        sasa_at = Dict{Atom,Float32}()

        for at in bb_atoms
            s = sasa(c.sasa_atoms, x -> x === at)
            sasa_at[at] = s
            bb_folded += s
        end
        for at in sc_atoms
            s = sasa(c.sasa_atoms, x -> x === at)
            sasa_at[at] = s
            sc_folded += s
        end

        bb_scale = iszero(bb_folded) ? 0.0f0 : bb_denatured / bb_folded
        sc_scale = iszero(sc_folded) ? 0.0f0 : sc_denatured / sc_folded

        bb_contrib = 0.0f0
        sc_contrib = 0.0f0

        for at in bb_atoms
            stype = record_surface_type(at)
            isnothing(stype) && continue
            σ = model_combination_rule(MTRecord, cos, stype)
            ΔASA = sasa_at[at] * (bb_scale - 1f0)
            bb_contrib += σ * ΔASA
        end

        for at in sc_atoms
            stype = record_surface_type(at)
            isnothing(stype) && continue
            σ = model_combination_rule(MTRecord, cos, stype)
            ΔASA = sasa_at[at] * (sc_scale - 1f0)
            sc_contrib += σ * ΔASA
        end

        residue_contributions_bb[i] = bb_contrib / 1000f0
        residue_contributions_sc[i] = sc_contrib / 1000f0
    end

    bb = sum(residue_contributions_bb)
    sc = sum(residue_contributions_sc)
    tot = bb + sc
    return MValue{MTRecord}(n, tot, bb, sc, residue_contributions_bb, residue_contributions_sc, cos)
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
    sasa_ats = sasa_particles(CreamerUnitedAtomRadii, atoms; unitcell)
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
    transfer_free_energy(::Type{MTRecord}, sasa_ats::SASA{CreamerUnitedAtomRadii}, cosolvent::AbstractString; kargs...)

Compute fixed-state transfer free energies using Record interaction potentials
and atom-resolved coarse-grained surface classes.

"""
function transfer_free_energy(
    ::Type{MTRecord},
    sasa_ats::SASA{CreamerUnitedAtomRadii},
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

            for at in bb_atoms
                stype = record_surface_type(at)
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

