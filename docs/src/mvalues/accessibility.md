```@meta
CollapsedDocStrings = true
```

!!! warning
    The implementation of this model is currently experimental. Parameterization details and 
    the interface might change.
    

# [The Accessibility model for transfer free energies](@id accessibility_model)

The `Accessibility` model revisits the Tanford additive transfer model from first principles, explicitly
accounting for the mutual shielding between the peptide backbone and the amino acid side chain. It requires
no experimental data beyond what is already used by the `AutonBolen` model, but redistributes the backbone
and side-chain contributions to transfer free energies (TFEs) and denaturation *m*-values in a way that is
consistent, simultaneously, with urea and with eight protecting osmolytes (TMAO, sarcosine, betaine, proline,
sorbitol, sucrose, glycerol, and trehalose).

**Motivation.** In the established (`AutonBolen`) model, the side-chain TFE of residue *aa* is obtained by
subtracting the full glycine TFE from the amino acid TFE:

```math
\text{TFE}^{sc}_{aa} = \text{TFE}_{aa} - \text{TFE}_\text{Gly}.
```

This removes exactly one backbone unit's worth of TFE, which is only correct if the backbone of *aa* is as
exposed as the backbone of glycine. The `MoeserHorinek` universal-backbone correction fixes this for the
backbone term alone, by referencing all backbone contributions to the glycine ASA, but it does not account
for the fact that the *side chain* is itself partly shielded by the backbone. As a result, it reproduces
urea denaturation well (once the glycine-activity correction is included) but systematically underestimates
the protective effect of protecting osmolytes, for which no analogous activity correction exists.

**Construction.** The `Accessibility` model decomposes the TFE of an amino acid into an *isolated* side
chain, an *isolated* backbone unit, and terminal (capping) groups, weighted by accessibility factors
``\alpha^{sc}_{aa}`` and ``\alpha^{bb}_{aa}`` (``\alpha=1``: fully exposed; ``\alpha=0``: fully shielded):

```math
\text{TFE}_{aa} = \alpha^{sc}_{aa}\, \text{TFE}^{i\text{-}sc}_{aa} + \alpha^{bb}_{aa}\, \text{TFE}^\text{bb} + \text{TFE}^\text{TG}.
```

Since glycine consists only of a backbone and terminal groups, ``\text{TFE}_\text{Gly} = \text{TFE}^\text{bb} + \text{TFE}^\text{TG}``,
which eliminates the (experimentally inaccessible) ``\text{TFE}^\text{TG}`` term and gives the TFE of the
*isolated* side chain:

```math
\text{TFE}^{i\text{-}sc}_{aa} = \frac{1}{\alpha^{sc}_{aa}} \left[ \text{TFE}_{aa} - \text{TFE}_\text{Gly} + \left(1 - \alpha^{bb}_{aa}\right)\text{TFE}^\text{bb} \right].
```

Using the apparent, activity-corrected transfer free energies of Auton & Bolen (``\text{TFE}^{sc,app}_{aa}``, with
the small amino-acid activity corrections ``\gamma_{aa}`` and ``\gamma_\text{Gly}``, the latter relevant only for
glycine in urea):

```math
\text{TFE}^{i\text{-}sc}_{aa} = \frac{1}{\alpha^{sc}_{aa}} \left[\text{TFE}^{sc,app}_{aa} + 
  \gamma_{aa} - \gamma_\text{Gly} + \left(1 - \alpha^{bb}_{aa}\right)\text{TFE}^\text{bb} \right].
```

The term ``\left(1-\alpha^{bb}_{aa}\right)\text{TFE}^\text{bb}`` corrects for the fact that subtracting
``\text{TFE}_\text{Gly}`` removes a *full* backbone unit, when only the exposed fraction of the backbone of
*aa* should be removed.

**Accessibility parameters.** ``\alpha^{sc}_{aa}`` and ``\alpha^{bb}_{aa}`` are obtained from accessible
surface area (ASA) calculations over a non-redundant protein structural database (CATH S20):

```math
\alpha^{sc}_{aa} = \frac{\text{ASA}^{sc}_{aa}}{\text{ASA}^{i\text{-}sc}_{aa}}, \qquad
\alpha^{bb}_{aa} = \frac{\text{ASA}^{bb}_{aa}}{\text{ASA}^\text{bb}}
```

where ``\text{ASA}^{sc}_{aa}`` and ``\text{ASA}^{bb}_{aa}`` are the side-chain and backbone ASAs of residue
*aa* in its normal (attached) context, and ``\text{ASA}^{i\text{-}sc}_{aa}`` and ``\text{ASA}^\text{bb}`` are
the ASAs of the side chain detached from the backbone and of the isolated backbone unit (equal to the
backbone ASA of glycine), respectively:

| Residue | ``\text{ASA}^{sc}`` (Å²) | ``\text{ASA}^{i\text{-}sc}`` (Å²) | ``\text{ASA}^{bb}_{aa}`` (Å²) | ``\text{ASA}^\text{bb}`` (Å²) | ``\alpha^{sc}_{aa}`` | ``\alpha^{bb}_{aa}`` |
|:--------|------:|------:|------:|------:|------:|------:|
| Ala | 70.88  | 135.35 | 49.01 | 88.89 | 0.5237 | 0.5514 |
| Phe | 187.59 | 248.85 | 38.66 | 87.37 | 0.7538 | 0.4425 |
| Leu | 159.14 | 219.49 | 37.32 | 88.26 | 0.7251 | 0.4229 |
| Ile | 159.82 | 220.11 | 35.32 | 87.57 | 0.7261 | 0.4034 |
| Val | 133.94 | 194.60 | 36.68 | 87.27 | 0.6883 | 0.4202 |
| Pro | 130.02 | 192.85 | 39.65 | 91.42 | 0.6742 | 0.4337 |
| Met | 161.46 | 222.97 | 39.91 | 87.86 | 0.7241 | 0.4543 |
| Trp | 230.51 | 291.63 | 37.45 | 88.01 | 0.7904 | 0.4255 |
| Gly |   0.00 |   0.00 | 87.48 | 87.48 | 1.0000 | 1.0000 |
| Ser |  85.27 | 148.71 | 45.90 | 87.71 | 0.5734 | 0.5234 |
| Thr | 117.74 | 179.11 | 39.44 | 87.03 | 0.6574 | 0.4532 |
| Tyr | 202.04 | 263.31 | 38.86 | 87.32 | 0.7673 | 0.4450 |
| Gln | 156.60 | 218.20 | 40.22 | 88.43 | 0.7177 | 0.4548 |
| Asn | 129.16 | 191.41 | 40.93 | 88.71 | 0.6748 | 0.4614 |
| Asp | 121.99 | 184.22 | 41.16 | 89.08 | 0.6622 | 0.4621 |
| Glu | 149.53 | 211.13 | 40.84 | 89.04 | 0.7082 | 0.4587 |
| His | 165.04 | 226.83 | 40.17 | 87.72 | 0.7276 | 0.4580 |
| Lys | 178.40 | 240.24 | 41.64 | 88.67 | 0.7426 | 0.4697 |
| Arg | 210.23 | 272.22 | 41.32 | 88.36 | 0.7723 | 0.4676 |
| Cys |  92.48 | 155.73 | 44.94 | 87.08 | 0.5939 | 0.5160 |

``\alpha^{bb}_{aa}`` is well below unity for all residues but glycine: the side chain shields the backbone
from direct solvent contact. Conversely, ``\alpha^{sc}_{aa}`` shows that the backbone shields at least
~20% of even the bulkiest side chains (Trp), an effect that previous models did not consider explicitly.

**Backbone accessibility is mechanism-dependent.** The geometric ratio ``\alpha^{bb}_{aa} = \text{ASA}^{bb}_{aa}/\text{ASA}^\text{bb}``
is appropriate when the cosolvent is excluded from, or interacts non-specifically with, the protein surface —
the case for protecting osmolytes. Urea, however, interacts with the backbone through hydrogen bonds, on the
face of the peptide unit opposite to the side chain; simulations show that the number of backbone–urea
hydrogen bonds is nearly independent of residue type, i.e., the side chain does *not* shield the backbone
from urea. The model therefore allows a mechanism parameter ``x`` per cosolvent:

```math
\alpha^{bb}_{aa}(x) = \left(\text{ASA}^{bb}_{aa}/\text{ASA}^\text{bb}\right)^{1-x}
```

with ``x=0`` recovering the geometric, ASA-based exposure (used for all protecting osmolytes) and ``x=1``
giving full backbone accessibility regardless of residue type, ``\alpha^{bb}_{aa}=1`` (used for urea). With
``x=1`` for urea, the backbone term of the model reduces exactly to the `MoeserHorinek` universal-backbone
term, recovering its predictive accuracy; with ``x=0``, side-chain shielding of the backbone is taken at
face value, which is appropriate for excluded-volume protectants.

**Final model.** Combining the isolated side-chain TFE with the universal backbone term gives the transfer
free energy of a protein:

```math
\text{TFE}^\text{prot} = \sum_{i=1}^{N_r} \left[
\left(\frac{\text{TFE}^{i\text{-}sc}_{aa}}{\text{ASA}^{i\text{-}sc}_{aa}}\right)\text{ASA}^{sc}_{aa} +
\left(\frac{\text{TFE}^\text{bb}}{\text{ASA}^\text{bb}_\text{Gly}}\right)\text{ASA}^{bb}_{aa}
\right]
```

where the sum runs over the ``N_r`` residues of the protein, and ``\text{ASA}^{sc}_{aa}`` and
``\text{ASA}^{bb}_{aa}`` are the actual (state-dependent) side-chain and backbone ASAs of each residue,
computed internally from the structure. The backbone term is identical to that of the `MoeserHorinek` model;
the side-chain term uses the isolated side-chain TFE and ASA, which differ from the quantities used by any
of the other models.

The model is available for all cosolvents in the Auton & Bolen parameterization (TMAO, Sarcosine, Betaine,
Proline, Sorbitol, Sucrose, Urea, Glycerol, and Trehalose):

```@example mvalue
using PDBTools
creamer_model = CreamerDenaturedModel(read_pdb(PDBTools.TESTPDB, "protein"))
m = mvalue(creamer_model, "urea"; model=Accessibility)
println("m-value: tot = $(m.tot), bb = $(m.bb), sc=$(m.sc)")
```

```@example mvalue
m = mvalue(creamer_model, "tmao"; model=Accessibility)
println("m-value: tot = $(m.tot), bb = $(m.bb), sc=$(m.sc)")
```

For protecting osmolytes, `Accessibility` predicts total *m*-values similar to `AutonBolen`, but with
side-chain contributions that are stabilizing and comparable to, or larger than, the backbone contribution —
in contrast to the backbone-dominated picture of `AutonBolen`. For urea, it reproduces the balanced
backbone/side-chain partition of `MoeserHorinek`. The predictions can differ from the models when
the ratio of backbone and side-chain accessibilities differ substantially from averages, or if the
exposed surface have very particular amino acid residue compositions. 