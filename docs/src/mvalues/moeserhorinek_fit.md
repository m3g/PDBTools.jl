```@meta
CollapsedDocStrings = true
```

# [The Gly-activity corrected Moeser&Horinek model](@id mh_gly)

The `MoeserHorinekFit` model is a variant of the Moeser & Horinek "universal backbone" transfer model in which a per-cosolvent
glycine activity coefficient correction is applied to the apparent side-chain transfer free energies (GTFEapp), extending the
physically motivated correction derived by Moeser & Horinek for urea to the other cosolvents available in the Auton & Bolen
parameterization.

**Background.** In the transfer model, side-chain transfer free energies (TFEs) are defined as the difference between the TFE
of an amino acid and that of glycine: ``\Delta\mu_\text{sc} = \Delta\mu_\text{aa} - \Delta\mu_\text{Gly}``. Because glycine
is highly soluble, its TFE carries a significant activity-coefficient correction that propagates to *all* 19 side-chain TFEs.
The commonly used GTFEapp values neglect this correction entirely; the GTFE\* values of Auton, Holthauzen & Bolen (2007)
include it, but with an algebraic error in the concentration-scale conversion that over-corrects by 25.8 cal mol⁻¹ M⁻¹.
Moeser & Horinek (2014) derived the corrected GTFE⁺ values (shift of −14.47 cal mol⁻¹ M⁻¹ relative to GTFEapp), which are
used in the `MoeserHorinek` model. For cosolvents other than urea, analogous experimental activity-coefficient data for
glycine are not generally available.

**The `MoeserHorinekFit` approach.** Rather than a first-principles correction, this model fits a single constant ``\gamma_G``
per cosolvent such that the resulting m-value predictions optimally reproduce those of the `AutonBolen` model. The corrected
side-chain TFE for residue *aa* in cosolvent *s* is:

```math
\Delta\mu_\text{sc,fit}(aa, s) = \Delta\mu_\text{sc,app}(aa, s) + \gamma_G(s)
```

where ``\Delta\mu_\text{sc,app}`` are the apparent GTFEapp values and ``\gamma_G`` is the fitted correction. The backbone is treated
as in the original Moeser & Horinek model (universal backbone: a single TFE per unit ASA for all residue types, using the backbone
ASA of Gly in isolated Gly-X-Gly sequences as reference). For urea, the fitted correction is ``\gamma_G \approx 14.1``
cal mol⁻¹ M⁻¹, consistent with the theoretical value of 14.47 cal mol⁻¹ M⁻¹ derived by Moeser & Horinek (2014), thus
validating the approach.

The `MoeserHorinekFit` model is available for all cosolvents parameterized in the `AutonBolen` model (TMAO, Sarcosine, Betaine,
Proline, Sorbitol, Sucrose, and Urea):

```@example mvalue
using PDBTools
creamer_model = CreamerDenaturedModel(read_pdb(PDBTools.TESTPDB, "protein"))
m = mvalue(creamer_model, "urea"; model=MoeserHorinekFit)
println("m-value: tot = $(m.tot), bb = $(m.bb), sc=$(m.sc)")
```

```@example mvalue
m = mvalue(creamer_model, "tmao"; model=MoeserHorinekFit)
println("m-value: tot = $(m.tot), bb = $(m.bb), sc=$(m.sc)")
```

The key distinction between the three models is how they handle the backbone and side-chain TFEs. In particular, the Moeser & Horinek models differ substantially in the prediction of the relative contribution of backbone and sidechains to the transfer free energies, relative to the (wrong) decomposition provided by the Auton and Bolen parameterization. 

| Model | Backbone | Side-chain TFEs | Cosolvents |
|:---|:---|:---|:---|
| `AutonBolen` | per-residue-type (established TM) | GTFE\* (activity-corrected, but erroneous) | all |
| `MoeserHorinek` | universal (Gly reference ASA) | GTFE⁺ (correctly corrected) | urea only |
| `MoeserHorinekFit` | universal (Gly reference ASA) | GTFEapp + fitted ``\gamma_G`` | all |
