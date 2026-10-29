module PDBTools

import Downloads
import JSON
import PrecompileTools
using Compat: @compat
using Format: format
using LinearAlgebra: norm, cross, dot
using Printf: @printf, @sprintf
using StaticArrays: SVector, SMatrix
using TestItems: @testitem, @testmodule
using OrderedCollections: OrderedDict
using InlineStrings: InlineString, String1, String15

#
# Core data structures and functions, defined in PDBToolsBase: the public names are
# re-exported (or declared public) below. Imported with `import`, such that PDBTools 
# can add methods to them (e.g. for selection strings).
#
import PDBToolsBase
import PDBToolsBase:
    # public names
    Atom, Chain, Element, MetaData, Model, NucleosideResidue, ProteinResidue, Residue, Segment, 
    Sequence, Structure, UnitCell,
    add_custom_field, add_element!, add_nucleoside_residue!, add_protein_residue!, 
    align, align!, alignment_movements, apply_alignment_transformation!,
    atomic_mass, atomic_number, atomic_symbol, beta, center_of_mass, chain, charge, closest,
    coordination_number, dihedral, dihedrals, distance, distances, 
    eachchain, eachmodel, eachresidue, eachsegment, 
    element, element_name, element_symbol, element_symbol_string, element_vdw_radius, elements,
    get_atoms, index, index_pdb, 
    isacidic, isaliphatic, isaromatic, isbackbone, isbasic, ischarged, ishydrophobic, isneutral,
    isnonpolar, isnucleoside, ispolar, isprotein, ispurine, ispyrimidine, issidechain, iswater,
    lattice_to_matrix, load, mass, matrix_to_lattice, maxmin, model, moveto!, name, 
    nucleoside_residues, occup, oneletter, pdb_element, position, positions, printatom, 
    protein_residues, read_mmcif, read_pdb, read_unitcell, read_xyz, remove_custom_elements!, 
    remove_custom_nucleoside_residues!, remove_custom_protein_residues!, residue, 
    residue_residue_distance, residuename, resname, resnum, rmsd, save, segname, set_position!, 
    threeletter, wrap, wrap_to_first, write_mmcif, write_pdb, write_xyz,
    # internal names
    AbstractStructuralElement, EachStructuralElement, StringType, backbone_atoms, not_side_chain_atoms,
    _resolve_unitcell

# Version
const VERSION = pkgversion(@__MODULE__)

# AtomsBase interface compatibility: deprecated until AtomsBase becomes stable
#import AtomsBase: atomic_number, atomic_symbol, atomic_mass, position
export atomic_number, atomic_symbol, atomic_mass, position, set_position!

export read_pdb, write_pdb, getseq, wget, edit!, oneletter, threeletter, residuename
export read_mmcif, write_mmcif
export read_xyz, write_xyz
export Atom, printatom, index, index_pdb, name, beta, occup, charge, pdb_element
export Structure, UnitCell
@compat public MetaData
@compat public test_file
@compat public Element, elements, ProteinResidue, protein_residues, NucleosideResidue, nucleoside_residues
@compat public align, align!, rmsd, alignment_movements, apply_alignment_transformation!
@compat public wrap, wrap_to_first, dihedrals, distances, coordination_number
export add_custom_field
export Residue, eachresidue, resname, residue, resnum, model, segname, zeta, zeta_check
export get_atoms
export Chain, eachchain, chain
export residue_ticks
export Segment, eachsegment
export Model, eachmodel
export positions, maxmin, distance, closest, dihedral, Ramachandran, set_phi!, set_psi!, extended_chain!, extended_chain
export element, mass, element_name, element_symbol, element_symbol_string, element_vdw_radius
export formula, stoichiometry
export Sequence
export isprotein, isnucleoside, isbackbone, issidechain, iswater
export isacidic, isaliphatic, isaromatic, isbasic, ischarged
export ishydrophobic, isneutral, isnonpolar, ispolar
export ispurine, ispyrimidine
export select_with_vmd
export visualize, open_browser
export VolumetricData, read_dx, write_dx

# Tools
export center_of_mass
export moveto!
export sasa_particles, sasa, SASA
export StandardAtomicRadii, CreamerUnitedAtomRadii, RichardsUnitedAtomRadii, CustomAtomicRadii
export read_unitcell, lattice_to_matrix, matrix_to_lattice
export hydrogen_bonds
export mvalue, CreamerDenaturedModel, save, load
@compat public parse_mvalue_server_sasa, gmx_delta_sasa_per_restype, delta_sasa_per_restype, creamer_delta_sasa

# Custom residue and element definitions
export custom_protein_residues!, add_protein_residue!, remove_custom_protein_residues!
export custom_nucleoside_residues!, add_nucleoside_residue!, remove_custom_nucleoside_residues!
export isnucleoside, ispurine, ispyrimidine
export custom_elements!, add_element!, remove_custom_elements!
export SIRAH

# Test files
const _extra_test_files = (
    "hbonds.pdb", "pbc.pdb", "no_pbc.pdb", "sirah.pdb",
    "1MJC_clean.pdb", "1MJC_native.pdb", "1MJC_straight.pdb",
    "2RN2_clean.pdb", "2RN2_native.pdb", "2RN2_straight.pdb", "2rmm_dimer.pdb",
)
"""
    PDBTools.test_file(name::AbstractString)

Returns the path of the test file `name`, used in the examples of the documentation 
(e.g. `read_pdb(PDBTools.test_file("structure.pdb"))`). The test files are downloaded 
on the first call (they are lazy artifacts, from https://github.com/m3g/MolSimTestFiles).

`test_file` is available after `using LazyArtifacts` (it is defined in a package extension).

"""
function test_file end

@testitem "test_file" begin
    using LazyArtifacts
    using PDBTools
    # files of PDBToolsBase and of PDBTools (lazy artifacts)
    for name in ("structure.pdb", "small.cif", "hbonds.pdb", "sirah.pdb", "1MJC_native.pdb")
        @test isfile(PDBTools.test_file(name))
    end
    @test all(isfile(PDBTools.test_file(name)) for name in PDBTools._extra_test_files)
end

# Analysis functions of residues
include("./struct_and_iterators/Residue.jl")

# Selection functions
include("./select/select.jl")
include("./select/select_with_vmd.jl")
include("./select/string_selections.jl")

#
# Input and output functions
#
include("./read_write/edit.jl")
include("./read_write/wget.jl")
include("./read_write/dx.jl")

include("./properties/formula.jl")
include("./properties/getseq.jl")

# Coordinates-dependent functions
include("./coordinates/contacts.jl")
include("./coordinates/dihedrals.jl")

# Miscellaneous tools
include("./tools/tools.jl")
include("./tools/sasa/sasa.jl")
include("./tools/sasa/cavity_exclusion.jl")
include("./tools/hydrogen_bonds.jl")
include("./tools/mvalue/mvalue.jl")
include("./tools/secondary_structure.jl")
include("./tools/visualize.jl")

# Custom element and residue definitions
include("./tools/custom_types.jl")

# Legacy compatibility

# Precompilation
include("./package_internals/precompile.jl")

end
