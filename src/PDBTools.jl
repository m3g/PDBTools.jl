module PDBTools

import Dates
import Downloads
import InteractiveUtils
import PrecompileTools
using Compat: @compat
using Format: format, generate_formatter
using LinearAlgebra: norm, cross, dot
using Printf: @printf, @sprintf
using StaticArrays: SVector, SMatrix
using TestItems: @testitem, @testmodule
using OrderedCollections: OrderedDict
using InlineStrings: InlineString, String1, String15
import MolSimToolkitShared: center_of_mass, wrap, dihedral

# Version
const VERSION = pkgversion(@__MODULE__)

# Default fixed String size for names, chains, segments, etc.
const StringType = String15

# AtomsBase interface compatibility: deprecated until AtomsBase becomes stable
#import AtomsBase: atomic_number, atomic_symbol, atomic_mass, position
export atomic_number, atomic_symbol, atomic_mass, position, set_position!

export read_pdb, write_pdb, getseq, wget, edit!, oneletter, threeletter, residuename
export read_mmcif, write_mmcif
export Atom, printatom, index, index_pdb, name, beta, occup, charge, pdb_element
export Structure
export add_custom_field
export Residue, eachresidue, resname, residue, resnum, model, segname, zeta, zeta_check
export get_atoms
export Chain, eachchain, chain
export residue_ticks
export Segment, eachsegment
export Model, eachmodel
export coor, maxmin, distance, closest, dihedral, Ramachandran
export element, mass, element_name, element_symbol, element_symbol_string, element_vdw_radius
export formula, stoichiometry
export Sequence
export select_with_vmd

# Tools
export center_of_mass
export moveto!
export sasa_particles, sasa
export read_unitcell, lattice_to_matrix, matrix_to_lattice
export hydrogen_bonds
export mvalue
@compat public parse_mvalue_server_sasa, gmx_delta_sasa_per_restype, delta_sasa_per_restype

# Custom residue and element definitions
export custom_protein_residues!, add_protein_residue!, remove_custom_protein_residues!
export custom_elements!, add_element!, remove_custom_elements!
export SIRAH

# Test files
const src_dir = @__DIR__
const test_dir = joinpath(src_dir, "../test")
const TESTPDB = joinpath(src_dir, "../test/structure.pdb")
const SMALLPDB = joinpath(src_dir, "../test/small.pdb")
const SIRAHPDB = joinpath(src_dir, "../test/sirah.pdb")
const TESTCIF = joinpath(src_dir, "../test/small.cif")
const CHAINSPDB = joinpath(src_dir, "../test/chains.pdb")
const BROKENCIF = joinpath(src_dir, "../test/broken.cif")
const LONG_CHAIN_STRING_CIF = joinpath(src_dir, "../test/long_chain_string.cif")
const LONG_NAME_STRING_CIF = joinpath(src_dir, "../test/long_name_string.cif")
const DIMERPDB = joinpath(src_dir, "../test/dimer.pdb")
const HETATMPDB = joinpath(src_dir, "../test/hetatm.pdb")
const TESTPBC = joinpath(src_dir, "../test/pbc.pdb")
const TESTNOPBC = joinpath(src_dir, "../test/no_pbc.pdb")

# Basic chemistry
include("./properties/elements.jl")
include("./properties/protein_residues.jl")

#
# Data structures and their iterators
#
include("./struct_and_iterators/Atom.jl")
include("./struct_and_iterators/Structure.jl")
include("./struct_and_iterators/iterators_base.jl")
include("./struct_and_iterators/Residue.jl")
include("./struct_and_iterators/Chain.jl")
include("./struct_and_iterators/Segment.jl")
include("./struct_and_iterators/Model.jl")

# Selection functions
include("./select/select.jl")
include("./select/select_with_vmd.jl")

#
# Input and output functions
#
include("./read_write/parsers.jl")
include("./read_write/read_pdb.jl")
include("./read_write/write_pdb.jl")
include("./read_write/read_mmcif.jl")
include("./read_write/write_mmcif.jl")
include("./read_write/edit.jl")
include("./read_write/wget.jl")

include("./properties/formula.jl")
include("./properties/getseq.jl")

# Coordinates-dependent functions
include("./coordinates/coor.jl")
include("./coordinates/distance.jl")
include("./coordinates/maxmin.jl")
include("./coordinates/contacts.jl")
include("./coordinates/dihedrals.jl")

# Miscellaneous tools
include("./tools/tools.jl")
include("./tools/sasa.jl")
include("./tools/hydrogen_bonds.jl")
include("./tools/read_unitcell.jl")
include("./tools/mvalue/mvalue.jl")
include("./tools/secondary_structure.jl")

# Custom element and residue definitions
include("./tools/custom_types.jl")

# Contact maps

# Legacy compatibility
include("./package_internals/legacy.jl")

# Precompilation
include("./package_internals/precompile.jl")

end
