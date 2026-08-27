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
parameterized directly on **atomic surface types**.

**Motivation.** Guinn et al. measured preferential interactions (osmometry and solubility) of urea
and glycine betaine (GB) with 44 small model compounds (amino acids, peptides, sugars, polyols,
aromatics, and salts), and fit the resulting chemical potential derivatives to a surface-additive
model over seven coarse-grained atom/surface types: aliphatic C, aromatic C, hydroxyl O, amide O,
carboxylate O, amide N, and cationic N. Because the surface types are transferable between model
compounds and proteins, the resulting interaction potentials can be applied directly to a protein
structure's atom-resolved solvent-accessible surface area (SASA), without ever forming an
amino-acid-level TFE.

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

Unlike the other transfer models, which pair a native structure with a *parametric* estimate of the
denatured-state backbone/side-chain ASA (the Creamer model, see
[`CreamerDenaturedModel`](protein_denaturation.md)), `MTRecord` is meant to be used with an explicit
*atomistic* model of the unfolded state, since it needs a real ``\Delta``ASA for every atom and
surface type, not just per-residue backbone/side-chain totals. The unfolded-state reference used
here is the fully extended (all-trans, ``\phi = \psi = 180°``) chain, obtained with
[`extended_chain`](@ref) - this is the same "extended β" reference state that Guinn et al.
themselves use in their own validation (Table S3 of the paper).

`MTRecordDenaturedModel` wraps a protein's native structure together with its extended chain (built
automatically) and the atom-resolved SASAs of both, computed once at construction:

```@docs
MTRecordDenaturedModel
mvalue(::MTRecordDenaturedModel, ::AbstractString)
```

```@example mvalue
model = MTRecordDenaturedModel(native_state)
m = mvalue(model, "urea")
```

The *m*-value is obtained directly from Eq. 4, using the real, atom-by-atom ``\Delta``ASA between
the extended and native chains (extended minus native). By default (`alpha=1.0`), the extended-chain
ASA is used as computed, with no adjustment; for urea, this reproduces the *m*-values reported in
Table S3 of Guinn et al. well. Glycine betaine has no equivalent extended-chain validation set in the
paper (Table S3 is urea-only), so predictions for `"betaine"` should be interpreted with more
caution; the `alpha` keyword is available to empirically rescale the extended-chain ASA, if needed
to better match experimental betaine data.
## Model details

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

**Interaction potentials.** The values of ``10^4\,\alpha_i`` (in `m⁻¹ Å⁻²`) currently implemented
are, for urea and betaine, from Table 1 of Guinn et al., and, for the other cosolvents, from the
surface-type potentials reported in the references below:

| Surface type | ``10^4\alpha_i`` (urea) | ``10^4\alpha_i`` (betaine) | ``10^4\alpha_i`` (TMAO) | ``10^4\alpha_i`` (proline) | ``10^4\alpha_i`` (trehalose) |
|:-------------|------:|------:|------:|------:|------:|
| Aliphatic carbon    | -1.1 |  3.0 | 17.8 |  5.3 | 22.4 |
| Aromatic carbon     | -8.9 | -23.0 | 22.6 | -9.2 |  5.9 |
| Hydroxyl oxygen     | -2.5 |  1.0 |  9.9 | -0.7 | -0.8 |
| Amide oxygen        | -8.7 | 28.0 |  9.3 | 14.5 | -19.6 |
| Carboxylate oxygen  | -4.0 | 29.0 | 52.3 | 16.6 | -28.2 |
| Amide nitrogen      | -3.2 | -20.0 | 11.6 | -11.8 | -4.7 |
| Cationic nitrogen   |  1.8 | -12.0 | -4.6 | -12.6 | 12.9 |

`"urea"` and `"betaine"` come from Guinn et al. (PNAS 2011, Table 1); `"tmao"` from
[doi:10.1016/j.bpj.2016.09.035](https://doi.org/10.1016/j.bpj.2016.09.035); `"proline"` from
[doi:10.1021/bi400683y](https://doi.org/10.1021/bi400683y); and `"trehalose"` from
[doi:10.1016/j.bpj.2015.05.037](https://doi.org/10.1016/j.bpj.2015.05.037). These are currently
the five cosolvents for which a full set of surface-type potentials is available.
