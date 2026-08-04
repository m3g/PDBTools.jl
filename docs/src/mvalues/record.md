```@meta
CollapsedDocStrings = true
```

!!! warning
    The implementation of this model is currently experimental. Parameterization details and 
    the interface might change.


# [The Record surface-type transfer model](@id record_model)

The `MTRecord` model implements the atom-resolved transfer model of Guinn, Pegram, Capp, Pollock,
and Record ([PNAS 2011](https://doi.org/10.1073/pnas.1109372108)), which explains why urea
denatures proteins while glycine betaine (a protecting osmolyte) stabilizes them. Unlike the
`AutonBolen`, `MoeserHorinek`, and `Accessibility` models, which are built from amino-acid transfer
free energies (TFEs) split into backbone and side-chain contributions, the `MTRecord` model is
parameterized directly on **atomic surface types**: it does not use amino-acid TFE data at all.

**Motivation.** Guinn et al. measured preferential interactions (osmometry and solubility) of urea
and glycine betaine (GB) with 44 small model compounds (amino acids, peptides, sugars, polyols,
aromatics, and salts), and fit the resulting chemical potential derivatives to a surface-additive
model over seven coarse-grained atom/surface types: aliphatic C, aromatic C, hydroxyl O, amide O,
carboxylate O, amide N, and cationic N. Because the surface types are transferable between model
compounds and proteins, the resulting interaction potentials can be applied directly to a protein
structure's atom-resolved solvent-accessible surface area (SASA), without ever forming an
amino-acid-level TFE.

**The additive model.** Guinn et al.'s Eq. 4 expresses the chemical potential derivative
``\mu_{23}/RT`` (essentially the preferential-interaction *m*-value per unit of exposed area) of a
compound as a sum over surface types ``i``, weighted by the corresponding ASA, plus an optional
ionic term:

```math
\frac{\mu_{23}}{RT} = \sum_i \alpha_i\, \text{ASA}_i + \sum_\text{ion} \nu_\text{ion}\, \beta_\text{ion}.
```

`PDBTools` uses only the surface-type sum (proteins are treated as neutral in this model). For a
surface type ``i`` and cosolvent, `PDBTools.model_combination_rule` returns the potential in
`cal mol⁻¹ Å⁻²`:

```math
\sigma_i = RT\,\alpha_i,
```

so that the contribution of a given atom to the transfer free energy (or *m*-value) is
``\sigma_i \times \text{ASA}`` (or ``\sigma_i \times \Delta\text{ASA}`` when comparing two states).

**Surface-type assignment.** `PDBTools.record_surface_type` maps each protein atom to one of the
seven surface types based on element and residue/atom name, following the assignment used by Guinn
et al.:

| Surface type | Assignment |
|:-------------|:-----------|
| Aliphatic carbon | Any carbon not classified as aromatic |
| Aromatic carbon | Ring carbons of Phe, Tyr, Trp, His |
| Hydroxyl oxygen | Ser OG, Thr OG1, Tyr OH |
| Amide oxygen | Backbone carbonyl O, Asn OD1, Gln OE1 |
| Carboxylate oxygen | Asp/Glu carboxylate oxygens, C-terminal OXT/OT1/OT2 |
| Amide nitrogen | Backbone amide N, Asn ND2, Gln NE2, and Trp NE1 (explicit exception) |
| Cationic nitrogen | Lys NZ, Arg NE/NH1/NH2, and His ND1/NE2 (explicit exception) |

Sulfur (Met, Cys) and hydrogen atoms are ignored, since their contribution to the ``\Delta``ASA of
folding is negligible in the reference dataset. These surface types are also exposed as selection
macro keywords, so that, for example, `select(atoms, "record_aromatic_carbon")` selects the atoms
assigned to the aromatic-carbon surface type.

**Interaction potentials.** The values of ``10^4\,\alpha_i`` (in `m⁻¹ Å⁻²`) currently implemented,
from Table 1 of Guinn et al., are:

| Surface type | ``10^4\alpha_i`` (urea) | ``10^4\alpha_i`` (betaine) |
|:-------------|------:|------:|
| Aliphatic carbon    | -1.1 |  3.0 |
| Aromatic carbon     | -8.9 | -23.0 |
| Hydroxyl oxygen     | -2.5 |  1.0 |
| Amide oxygen        | -8.7 | 28.0 |
| Carboxylate oxygen  | -4.0 | 29.0 |
| Amide nitrogen      | -3.2 | -20.0 |
| Cationic nitrogen   |  1.8 | -12.0 |

Only `"urea"` and `"betaine"` are currently supported, since these are the two cosolvents for which
Guinn et al. report a full set of surface-type potentials.

## Fixed-state transfer free energies

The transfer free energy of a single (fixed) structure can be computed directly from its
atom-resolved SASA with `transfer_free_energy`, by requesting the `MTRecord` model:

```@example mvalue
using PDBTools
native_state = read_pdb(PDBTools.MJC_NATIVE, "protein")
tfe = transfer_free_energy(native_state, "urea"; model=MTRecord)
```

As with the other models, the result is split into backbone and side-chain contributions, here
obtained by summing the atom-resolved surface-type contributions of backbone and side-chain atoms,
respectively.

## Denaturation *m*-values

Because `MTRecord` requires atom-resolved SASAs on both the folded and the denatured state, and the
denatured state has no explicit atomic structure, *m*-values of denaturation are computed by
wrapping the atoms (or an existing `CreamerDenaturedModel`) in a `MTRecordDenaturedModel`, which
reuses the Creamer backbone/side-chain denatured-state ASA estimates (as in
[`CreamerDenaturedModel`](protein_denaturation.md)) to scale the folded, atom-resolved
``\Delta``ASA of each residue:

```@docs
MTRecordDenaturedModel
mvalue(::MTRecordDenaturedModel, ::AbstractString)
```

```@example mvalue
model = MTRecordDenaturedModel(native_state)
m = mvalue(model, "urea")
println("m-value: tot = $(m.tot), bb = $(m.bb), sc = $(m.sc)")
```

The `alpha` keyword of `mvalue` (default `1.15`) rescales the Creamer denatured-state backbone and
side-chain ASA estimates before they are combined with the Record surface-type potentials. This
empirical correction accounts for the fact that Creamer's denatured-state ASAs were calibrated for
use with amino-acid-based additive models (`AutonBolen`), and are not, without adjustment, on the
same scale as the atom-resolved ``\Delta``ASA required by `MTRecord`. The default value was found to
reproduce the urea benchmarks of Guinn et al. (Table S3) well; for other cosolvents or particularly
unusual proteins, `alpha` (and the denaturation `type`, `1`, `2`, or `3`, of `MTRecordDenaturedModel`)
may need to be adjusted to match experimental data, as illustrated in the package tests for betaine.
