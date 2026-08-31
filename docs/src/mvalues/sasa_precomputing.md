```@meta
CollapsedDocStrings = true
```

# High-throughput calculations with precomputed SASAs

The m-value and transfer free energy calculations both require computing the solvent
accessible surface area (SASA) of the protein atoms. Every model except `MTRecord` (see
[the Record model](@ref record_model)) is parameterized with Creamer united-atom radii;
`MTRecord` was parameterized against, and expects, Richards united-atom radii instead (to
match SurfaceRacer, which was used to compute its own ASA reference values). By default,
each call to `mvalue` or `transfer_free_energy` recomputes the SASA from scratch, using
whichever radii the requested `model` needs, which dominates the runtime. When only the
cosolvent, selection, or model changes between calls — but the atomic coordinates stay
fixed — the SASA can be computed once and reused.

## Computing the SASA with Creamer radii

Pass `CreamerUnitedAtomRadii` as the first argument to `sasa_particles` to obtain a
`SASA{CreamerUnitedAtomRadii}` object that is compatible with both `mvalue` and
`transfer_free_energy`:

```@example mvalue
using PDBTools
native_state = read_pdb(PDBTools.MJC_NATIVE, "protein")
desnat_state = read_pdb(PDBTools.MJC_DESNAT, "protein")
sasa_native  = sasa_particles(CreamerUnitedAtomRadii, native_state)
```

```@example mvalue
sasa_desnat  = sasa_particles(CreamerUnitedAtomRadii, desnat_state)
```

## Reusing precomputed SASAs for m-values

Pass the two `SASA{CreamerUnitedAtomRadii}` objects directly to `mvalue` instead of the
atom arrays. All keyword arguments (`model`, `sel`, `backbone`, `sidechain`, `parallel`)
work identically:

```@example mvalue
m_urea = mvalue(sasa_native, sasa_desnat, "urea"; model=MoeserHorinek)
```

The same precomputed SASAs can be queried immediately without repeating the geometry computation,
with a different cosolvent
```@example mvalue
m_tmao   = mvalue(sasa_native, sasa_desnat, "tmao")
```

or a residue selection,

```@example mvalue
m_acidic = mvalue(sasa_native, sasa_desnat, "urea"; sel="acidic")
```

## Reusing precomputed SASAs for transfer free energies

For `transfer_free_energy`, a single structure is used, so only one SASA is precomputed:

```@example mvalue
sasa_native2 = sasa_particles(CreamerUnitedAtomRadii, native_state)
```

```@example mvalue
t_urea = transfer_free_energy(sasa_native2, "urea")
```

```@example mvalue
t_tmao = transfer_free_energy(sasa_native2, "tmao")
```

## Reusing precomputed SASAs with the Record model (`MTRecord`)

`MTRecord` needs Richards, not Creamer, united-atom radii, so precompute the SASA with
`RichardsUnitedAtomRadii` instead:

```@example mvalue
sasa_native_richards = sasa_particles(RichardsUnitedAtomRadii, native_state)
```

The resulting `SASA{RichardsUnitedAtomRadii}` object can be reused across cosolvents, exactly
as with the Creamer-radii SASAs above, by passing `model=MTRecord`:

```@example mvalue
t_urea_record = transfer_free_energy(sasa_native_richards, "urea"; model=MTRecord)
```

```@example mvalue
t_betaine_record = transfer_free_energy(sasa_native_richards, "betaine"; model=MTRecord)
```

The same applies to `MTRecordDenaturedModel`, which is itself just a pair of
`SASA{RichardsUnitedAtomRadii}` objects (native and fully-extended chain) computed once at
construction and reused for every cosolvent passed to `mvalue` — see
[Denaturation *m*-values](@ref record_model).

## Enforced radii compatibility

Passing a `SASA` object computed with the wrong radii parameterization for the requested model
raises an informative error (naming the model and the radii type it expects), ensuring that
m-values and transfer free energies are never silently computed from incompatible surface areas.
For example, `MTRecord` requires `RichardsUnitedAtomRadii`, while every other model requires
`CreamerUnitedAtomRadii`; passing one where the other is expected raises an `ArgumentError`
instead of silently returning a wrong result.
