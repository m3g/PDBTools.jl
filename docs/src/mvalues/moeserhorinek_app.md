```@meta
CollapsedDocStrings = true
```

!!! warn
    This is experimental and might be discontinued.

# [The Universal Backbone Moeser & Horinek model](@id mh_app)

The `MoeserHorinekApp` model applies the Moeser & Horinek universal backbone to the Auton & Bolen parameterization,
recalculating the side-chain transfer free energies (SC TFEs) consistently with the universal backbone. No
glycine-activity correction is applied.

**Background.** In the Auton&Bolen transfer model (TM), the side-chain TFE of residue *aa* is obtained as:

```math
\Delta G^\text{est}_\text{sc}(aa) = \text{TFE}_{aa} - \text{TFE}_\text{Gly}
```

Subtracting glycine's full TFE removes exactly one backbone unit's worth of TFE, because glycine's isolated ASA
is 88.1 Å² and the backbone contribution per unit ASA in the universal backbone is
``\text{TFE}_\text{bb} / 88.1``. For any other residue, however, the backbone ASA ``A^\text{bb}_{aa} \neq 88.1``
Å², so the backbone contribution is ``A^\text{bb}_{aa} \times \text{TFE}_\text{bb} / 88.1``, which is
**not** equal to ``\text{TFE}_\text{Gly}`` (unless the residue is glycine itself). The apparent GTFEapp values
thus embed an physically unjustified backbone subtraction for every non-glycine residue.

**Correction.** The proper universal-backbone SC TFE proposed by Moeser&Horinek, which leads to a universal backbone
contribution for all residues is:

```math
\Delta G^\text{UB}_\text{sc}(aa) = \text{TFE}_{aa} - A^\text{bb}_{aa} \times \frac{\text{TFE}_\text{bb}}{88.1}
= \Delta G^\text{est}_\text{sc}(aa) + \text{TFE}_\text{bb} \left(1 - \frac{A^\text{bb}_{aa}}{88.1}\right)
```

This correction is amino-acid-specific. For example (urea, ``\text{TFE}_\text{bb} = -39`` cal mol⁻¹ M⁻¹):

| Residue | ``A^\text{bb}`` (Å²) | Correction (cal mol⁻¹ M⁻¹) |
|:--------|---------------------:|-----------------------------:|
| Gly | 88.1 | 0 |
| Ala | 46.2 | +19.1 |
| Ile | 30.9 | +25.3 |

The `MoeserHorinekApp` model applies this residue-specific correction to the GTFEapp values from Auton & Bolen,
uses the universal backbone (``A^\text{bb}_{aa} \times \text{TFE}_\text{bb} / 88.1``), and *does not* apply
any glycine-activity correction. As a result, for urea specifically, `MoeserHorinekApp` is less accurate than
`MoeserHorinek`, which includes the Gly-activity correction derived by Moeser & Horinek (2014).

The model is available for all cosolvents in the Auton & Bolen parameterization (TMAO, Sarcosine, Betaine,
Proline, Sorbitol, Sucrose, Urea, Glycerol, and Trehalose):

```@example mvalue
using PDBTools
creamer_model = CreamerDenaturedModel(read_pdb(PDBTools.TESTPDB, "protein"))
m = mvalue(creamer_model, "urea"; model=MoeserHorinekApp)
println("m-value: tot = $(m.tot), bb = $(m.bb), sc=$(m.sc)")
```

```@example mvalue
m = mvalue(creamer_model, "tmao"; model=MoeserHorinekApp)
println("m-value: tot = $(m.tot), bb = $(m.bb), sc=$(m.sc)")
```

The key distinction between the three models is how they handle the backbone and side-chain TFEs:

| Model | Backbone | Side-chain TFEs | Gly-activity correction | Cosolvents |
|:---|:---|:---|:---|:---|
| `AutonBolen` | per-residue-type (established TM) | GTFEapp | none | all |
| `MoeserHorinek` | universal (Gly reference ASA) | GTFE⁺ (correctly corrected for urea) | yes | urea only |
| `MoeserHorinekApp` | universal (Gly reference ASA) | GTFEapp recalculated via Eq. above | none | all |
