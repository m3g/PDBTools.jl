```@meta
CollapsedDocStrings = true
```

# [The Universal Backbone Moeser & Horinek model](@id mh_app)

The `MoeserHorinekApp` model applies the Moeser & Horinek universal backbone to the Auton & Bolen parameterization,
recalculating the side-chain transfer free energies (SC TFEs) consistently with the universal backbone. No
glycine-activity correction is applied.

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

## Background

In the Auton&Bolen transfer model (TM), the contribution to the protein transfer free energy of
an amino acid residue $$r$$ undergoing a change is surface area exposure to the solvent is
```math
\Delta TFE(r) = \Delta ASA^{sc}(r) \left(\frac{TFE^{sc}_{aa}}{iASA^{sc}_{aa}}\right) +
         \Delta ASA^{bb}(r) \left(\frac{TFE_{Gly}}{iASA^{bb}_{aa}}\right)
```
where $$\Delta ASA^{sc}(r)$$ is the change in surface area exposure of the side chain,
$$\Delta ASA^{bb}(r)$$ that of the backbone. $$TFE^{sc}_{aa}$$ is the transfer free energy of the 
amino acid type $$aa$$, and $$iASA^{sc}_{aa}$$ stands for the surface area os the isolated
residue. The crucial part is that for the backbone, Auton and Bolen use the transfer
free energy of glycine, $$TFE_{Gly}$$ divided by the isolated-residue surface area of
the backbone *of the amino acid type*, $$iASA^{bb}_{aa}$$. Given that the isolated 
backbone surface are is smaller than that of glycine (because of side-chain protection), 
the term $$TFE_{Gly}/iASA^{bb}_{aa}$$
overestimates the contribution of the backbone. Also, by using a amino acid type specific
transfer free energy per unit area of the backbone, the model breaks strict group
additivity, as the backbone of all residues are chemically identical.     

Moeser and Horinek noticed this inconsitency and proposed that the proper, and simpler 
evaluation of an amino acid contribution should be
```math
\Delta TFE(r) = 
         \Delta ASA^{sc}(r) \left(\frac{TFE^{sc}_{aa}}{iASA^{sc}_{aa}}\right) +
         \Delta ASA^{bb}(r) \left(\frac{TFE_{Gly}}{iASA_{Gly}}\right)
```
where the crucial difference is the use of the full isolated-Glycine ASA in the denominator
of the last term, with two consenquences: 1) The contributions of the backbones are 
relatively reduced, and 2) all backbone contributions are defined by the same chemical
nature restoring group additivity. They named this model the "Universal Backbone" 
additive model.

The two models can provide the same estimates for the total transfer free energy of an amino acid, 
but this is dependent then on how the side chain contributions are defined from experimental data.
For consistency with the Auton and Bolen model, the apparent side chain transfer free energies must 
be
```math
TFE^{sc-AB}_{aa} = TFE_{aa} - TFE_{Gly}
```
where $$aa$$ or $$Gly$$ subscript usually refer to tri-peptide Gly-X-Gly residues. With this definiton,
a fully-exposed backbone is subtracted from the amino acid $$TFE$$, ignoring the fact that the
side chain itself protects the backbone. The Auton and Bolen model correct this by later defining
residue type-specific backbone contributions, as shown above.  

The alternative model is to incorporate directly into $$TFE^{sc}_{aa}$$ the fact that the side chain
protects the backbone differently for each amino acid type. That is, 
```math
TFE^{sc-MH}_{aa} = TFE_{aa} - \left(\frac{iASA^{aa}_{bb}}{iASA^{Gly}}\right) TFE_{Gly}
```
where the term in parenthesis is the fraction of the backbone that is exposed to solvent in 
the isolated amino acid type $$aa$$. 

By replacing $$TFE_{aa}$$ from the equation of $$TFE^{sc-AB}_{aa}$$ into the equation for $$TFE^{sc-MH}_{aa}$$ 
we get
```math
TFE^{sc-MH}_{aa} =TFE^{sc-AB}_{aa} + \left[ 1 - \left(\frac{iASA^{aa}_{bb}}{iASA^{Gly}}\right)\right] TFE_{Gly}
```
which is a side-chain transfer free energy which takes into account that part of the backbone
is protected from the solvent, differently, for each amino acid type. This allows computing universal backbone
model (Moeser & Horinek) side-chain contributions from the Auton & Bolen tabulated side-chain contributions
and isolated amino acid and backbone areas, and the fact that $$TFE_{Gly}$$ is constant for each cosolvent.

This correction is amino-acid-specific. For example (urea, ``\text{TFE}_\text{bb} = -39`` cal mol⁻¹ M⁻¹):

| Residue | ``A^\text{bb}_{aa}`` (Å²) | Correction (cal mol⁻¹ M⁻¹) |
|:--------|---------------------:|-----------------------------:|
| Gly | 88.1 | 0 |
| Ala | 46.2 | +19.1 |
| Ile | 30.9 | +25.3 |

The `MoeserHorinekApp` model applies this residue-specific correction to the GTFEapp values from Auton & Bolen,
uses the universal backbone (``A^\text{bb}_{aa} \times \text{TFE}_\text{Gly} / 88.1``), and *does not* apply
any glycine-activity correction. As a result, for urea specifically, `MoeserHorinekApp` is less accurate than
`MoeserHorinek`, which includes the Gly-activity correction derived by Moeser & Horinek (2014).
