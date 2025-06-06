PDBTools.jl Changelog
===========================

  
[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-experimental]: https://img.shields.io/badge/Experimental-yellow.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-fix]: https://img.shields.io/badge/Fix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

Version 3.1.1
-------------
- ![ENHANCEMENT][badge-enhancement] better handling of alternate conformations of protein residues. 
- ![INFO][badge-info] Set version to 3.1.1

Version 3.1.0
-------------
- ![FEATURE][badge-feature] Selection syntax now supports parenthesis and the shortcut for `or`, as in `resname ARG CYS ALA`. 
- ![INFO][badge-info] CYS is now classified as a polar residue, following the definition of VMD.
- ![INFO][badge-info] GLY is now classified as a polar residue, following the definition of VMD.
- ![INFO][badge-info] Set version to 3.1.0

Version 3.0.0
-------------
- ![FEATURE][badge-feature] `write_mmcif` and `write_pdb` support a third positional argument, `String` or `Function`, to write a selection of the atoms.
- ![ENHANCEMENT][badge-enhancement] `select_with_vmd` returns a vector of `Atom`s if the input was a vector of `Atom`s. 
- ![BREAKING][badge-breaking] `select_with_vmd`  returns a vector of `Atom`s if the input was a vector of `Atom`s. 
- ![BREAKING][badge-breaking] In all functions that accepted the `only` keyword parameter to define selections with a Julia function, the keyword parameter was dropped and now the function can be provided as a positional second argument, or third argument for `write_pdb` and `write_mmcif` functions.
- ![BREAKING][badge-breaking] `get_seq` functions do not support anymore the input of the file name. An array of `Atom`s must be provided. This is because internally it would be required to recognize PDB of mmCIF functions.
- ![INFO][badge-info] Dropped Julia 1.9 support (minimum requirement is 1.10)
- ![INFO][badge-info] Set version to 3.0.0 