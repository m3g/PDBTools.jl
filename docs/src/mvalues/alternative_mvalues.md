```@meta
CollapsedDocStrings = true
```

# Alternative m-value calculations

The following functions can be used to compute *m*-values from the variation of the SASA per residue 
type, which allow the use of external tools to compute the SASA. This is used mostly for testing 
purposes. The functions allow the use of SASAs obtained directly from the [Auton & Bolen server](http://best.bio.jhu.edu/mvalue/), or 
from Gromacs SASA calculations. The `creamer_delta_sasa` function uses the same random coil models of the 
server and same atom radii, providing the same results for unfolding m-values.

```@docs
PDBTools.mvalue_delta_sasa
PDBTools.creamer_delta_sasa
PDBTools.delta_sasa_per_restype
PDBTools.parse_mvalue_server_sasa
PDBTools.gmx_delta_sasa_per_restype
```
