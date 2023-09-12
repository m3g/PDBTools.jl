module PDBTools

using Parameters
using Printf
using InteractiveUtils
using Formatting
using StaticArrays
using LinearAlgebra: norm
using TestItems
import Dates

# AtomsBase interface compatibility
import AtomsBase: atomic_number, atomic_symbol, atomic_mass, position
export atomic_number, atomic_symbol, atomic_mass, position

export readPDB, writePDB, getseq, wget, edit!, oneletter, threeletter, residuename
export Atom, printatom, index, index_pdb, name, beta, occup
export Residue, eachresidue, resname, residue, resnum, chain, model, segname
export residue_ticks
export coor, maxmin, distance, closest
export element, mass, element_name
export formula, stoichiometry
export Sequence

const TESTPDB = "$(@__DIR__)/../test/structure.pdb"

# Basic chemistry
include("./elements.jl")
include("./protein_residues.jl")

#
# Data structures
#
include("./Atom.jl")
include("./Residue.jl")

include("./Indexes_mmCIF_fields.jl")

# Selection functions
include("./select.jl")

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
