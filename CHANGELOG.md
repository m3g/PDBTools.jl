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

Version 3.19.0
--------------
- ![FEATURE][badge-feature] Add `parallel` option to contact and distance map calcualtions, set to false by default.
- ![ENHANCEMENT][badge-enhancement] Implement fast contact and distance map calculations using cell lists. 
- ![INFO][badge-info] Document properly all property getters. 

Version 3.18.0
--------------
- ![FEATURE][badge-feature] Numeric selection keywords accept "to" to define ranges, as in `residue 1 to 5`. 

Version 3.17.0
--------------
- ![FEATURE][badge-feature] Implement `transfer_free_energy` function to compute the TFE of a single structure to different solutions. 
- ![ENHANCEMENT][badge-enhancement] MValue object carries the `cosolvent` field, and shows the cosolvent.

Version 3.16.1
--------------
- ![ENHANCEMENT][badge-enhancement] Better support for mmCIF files when `ATOM` strings are found outside atom descriptors. 

Version 3.16.0
--------------
- ![FEATURE][badge-feature] Compute m-value of protein denaturation using Creamer denatured models.
- ![FEATURE][badge-experimental] Introduce `Structure` object, to contain atom arrays and meta-data.
- ![FEATURE][badge-experimental] `Atom` object equality is supported (all fields must be identical).

Version 3.15.1
--------------
- ![INFO][badge-info] Use `String15` to store names, chains, etc, to support very large CIF files.
- ![INFO][badge-info] The maximum string length for names, chains, etc, is defined by the constant global StringType.

Version 3.15.0
--------------
- ![FEATURE][badge-feature] Allow updating element properties when adding custom element.
- ![FEATURE][badge-feature] `sasa_particles(SIRAH, args...; kargs...)` computes SASAs with SIRAH parameters.
- ![ENHANCEMENT][badge-enhancement] Update vdW radii and masses of SIRAH elements.
- ![ENHANCEMENT][badge-enhancement] Add c-terminal and n-terminal SIRAH residues to custom residue list.

Version 3.14.0
--------------
- ![FEATURE][badge-feature] Append SIRAH backbone atoms to allow selection of backbone and sidechains with SIRAH pdbs.
- ![INFO][badge-info] Use julia-cations/cache

Version 3.13.0
--------------
- ![FEATURE][badge-feature] `ContactMap` objects can be summed or subtracted, to obtain the effect of transformations on contact maps.

Version 3.12.0
--------------
- ![FEATURE][badge-feature] add `parallel` keyword argument to `mvalue` function to enable parallel computation of residue contributions.
- ![ENHANCEMENT][badge-enhancement] Improve performance of mvalue calculations by avoiding unnecessary allocations.

Version 3.11.3
--------------
- ![ENHANCEMENT][badge-enhancement] Small improvements in performance, particularly of the mvalue calculations. 
- ![ENHANCEMENT][badge-enhancement] Cosolvent selection in `mvalue` is now case-insensitive.

Version 3.11.2
--------------
- ![ENHANCEMENT][badge-enhancement] The `MValue` data structure now contains a type parameter that indicates the model used, and a `Int` field with the number of residues in the selection.

Version 3.11.1
--------------
- ![ENHANCEMENT][badge-enhancement] Use `InlineStrings.String7` for `chain` in `Atom` struct to support longer chain IDs from mmCIF files.
- ![ENHANCEMENT][badge-enhancement] Faster `in` implementation for containment of atoms in residues, chains, models, segments.

Version 3.11.0
--------------
- ![FEATURE][badge-feature] add `stride_run` and `dssp_run` functions to compute secondary structure from vectors of atoms.

Version 3.10.0
--------------
- ![FEATURE][badge-feature] add `mvalue` function for computing *m*-values from protein structures.
- ![FEATURE][badge-feature] `in` is part of supported interface for iterators (e. g. `atoms[1] in residues[1]`)

Version 3.9.0
-------------
- ![FEATURE][badge-feature] add `get_atoms` getter for fetching the vector of atoms of residues, models, chains, segments.
- ![FEATURE][badge-feature] `hydrogen_bonds` accepts multiple pairs of selections, and is made stable.
- ![INFO][badge-info] Documentation updates.

Version 3.8.0
-------------
- ![FEATURE][badge-experimental] `hydrogen_bonds` function to compute h-bonds in structures (with Hydrogen atoms added).
- ![INFO][badge-info] Add tests to computing contact maps with PBCs.

Version 3.7.0
-------------
- ![FEATURE][badge-feature] `set_position!(::Atom, ::Union{Tuple,AbstractVector})` to conveniently set coordinates from a vector or tuple of coordinates.
- ![FEATURE][badge-feature] Support for periodic boundary conditions in `sasa_particles` with `unitcell` keyword argument.
- ![FEATURE][badge-feature] Add `read_unitcell`, `lattice_to_matrix`, and `lattice_to_matrix` functions to read and convert unitcells written in PDB or mmCIF files.

Version 3.6.0
-------------
- ![FEATURE][badge-feature] `output_dots` option to `sasa_particles` to return the dots of the surface. 
- ![ENHANCEMENT][badge-enhancement] return more compact SASA object. Rename `atomic_sasa` to `sasa_particles` (keep alias for compatibility).
- ![ENHANCEMENT][badge-enhancement] make the `AtomDots` of sasa calculation a continguous memory block, to reduce GC pressure.
- ![INFO][badge-info] Set default `n_dots` of `atomic_sasa` computations to 512 instead of 500 such that it is a multiple of 16, the default `N_SIMD`.
- ![INFO][badge-info] add doc examples.

Version 3.5.5
-------------
- ![ENHANCEMENT][badge-enhancement] use Fibonacci lattice to generate exactly `n_dots` the dots in the sphere. Use default `n_dots=500`.

Version 3.5.4
-------------
- ![INFO][badge-info] Provide a better error message when a field in a CIF file does not fit in the fixed-size string to which it must be assigned.

Version 3.5.3
-------------
- ![ENHANCEMENT][badge-enhancement] precompile SASA functions. 

Version 3.5.2
-------------
- ![ENHANCEMENT][badge-enhancement] faster SASA calculation.

Version 3.5.1
-------------
- ![ENHANCEMENT][badge-enhancement] improve selection syntax for non-contiguous indexing in `sasa`. 

Version 3.5.0
-------------
- ![FEATURE][badge-feature] `atomic_sasa` and `sasa` functions to compute solvent accessible surface area.
- ![FEATURE][badge-feature] Add `vdw_radius` as element property, and `element_vdw_radius` function to fetch them.
- ![INFO][badge-info] The element symbol is returned as a `InlineStrings.String3`. 
- ![INFO][badge-info] Masses are returned as `Float32` values instead of `Float64`.

Version 3.4.0
-------------
- ![FEATURE][badge-feature] `Plots.scatter(::Ramachandran)` to plot Ramachandran objects.
- ![FEATURE][badge-feature] `Ramachandran` method and struct, to compute Ramachandran plots.
- ![FEATURE][badge-feature] `dihedral` method to compute directly dihedral given four Atom objects.

Version 3.3.0
-------------
- ![FEATURE][badge-feature] `zeta` and `zeta_check` functions to check for Calpha chirality of protein residues. 

Version 3.2.0
-------------
- ![FEATURE][badge-feature] Selection by coordinate values (e. g. `select(atoms,"x > 0")`) 

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