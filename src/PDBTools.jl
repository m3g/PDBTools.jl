module PDBTools

import Dates
import Downloads
import InteractiveUtils
import PrecompileTools
using Format: format
using LinearAlgebra: norm
using Parameters: @unpack
using Printf: @printf, @sprintf
using StaticArrays: SVector
using TestItems: @testitem
using OrderedCollections: OrderedDict
using InlineStrings: String3, String7

# Version
const VERSION = pkgversion(@__MODULE__)

# AtomsBase interface compatibility: deprecated until AtomsBase becomes stable
#import AtomsBase: atomic_number, atomic_symbol, atomic_mass, position
export atomic_number, atomic_symbol, atomic_mass, position

export read_pdb, write_pdb, getseq, wget, edit!, oneletter, threeletter, residuename
export read_mmcif, write_mmcif
export Atom, printatom, index, index_pdb, name, beta, occup, charge, pdb_element
export add_custom_field
export Residue, eachresidue, resname, residue, resnum, chain, model, segname
export residue_ticks
export coor, maxmin, distance, closest
export element, mass, element_name, element_symbol, element_symbol_string
export formula, stoichiometry
export Sequence
export select_with_vmd

# Tools
export center_of_mass
export moveto!

# Custom residue and element definitions
export custom_protein_residues!, add_protein_residue!, remove_custom_protein_residues!
export custom_elements!, add_element!, remove_custom_elements!
export SIRAH 

# Test files
const TESTPDB = joinpath(@__DIR__,"../test/structure.pdb")
const SMALLPDB = joinpath(@__DIR__,"../test/small.pdb")
const SIRAHPDB = joinpath(@__DIR__,"../test/sirah.pdb")
const TESTCIF = joinpath(@__DIR__,"../test/small.cif")

# Basic chemistry
include("./elements.jl")
include("./protein_residues.jl")

#
# Data structures
#
include("./Atom.jl")
include("./Residue.jl")

# Selection functions
include("./select.jl")
include("./select_with_vmd.jl")

#
# Input and output functions
#
include("./parsers.jl")
include("./formula.jl")
include("./read_pdb.jl")
include("./write_pdb.jl")
include("./read_mmcif.jl")
include("./write_mmcif.jl")
include("./edit.jl")
include("./getseq.jl")
include("./coor.jl")
include("./distance.jl")
include("./maxmin.jl")
include("./wget.jl")
include("./tools.jl")

# Custom element and residue definitions
include("./custom_types.jl")

# Legacy compatibility
include("./legacy.jl")

# Precompilation
include("./precompile.jl")

end
