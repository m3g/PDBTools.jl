module PDBTools

import Dates
import Downloads
import InteractiveUtils
using Formatting: format
using LinearAlgebra: norm
using Parameters: @unpack
using Printf: @printf, @sprintf
using StaticArrays: SVector
using TestItems: @testitem

# AtomsBase interface compatibility
import AtomsBase: atomic_number, atomic_symbol, atomic_mass, position
export atomic_number, atomic_symbol, atomic_mass, position

export readPDB, writePDB, getseq, wget, edit!, oneletter, threeletter, residuename
export Atom, printatom, index, index_pdb, name, beta, occup, custom_field
export Residue, eachresidue, resname, residue, resnum, chain, model, segname
export residue_ticks
export coor, maxmin, distance, closest
export element, mass, element_name, element_symbol
export formula, stoichiometry
export Sequence
export select_with_vmd

const TESTPDB = "$(@__DIR__)/../test/structure.pdb"
const SMALLPDB = "$(@__DIR__)/../test/small.pdb"

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
include("./formula.jl")
include("./check_mmCIF.jl")
include("./read_atom.jl")
include("./write_atom.jl")
include("./readPDB.jl")
include("./edit.jl")
include("./writePDB.jl")
include("./getseq.jl")
include("./coor.jl")
include("./distance.jl")
include("./maxmin.jl")
include("./wget.jl")

end
