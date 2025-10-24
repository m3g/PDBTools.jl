# PDBTools.jl

A lightweight and flexible Julia package for handling protein structure files with ease.

## Installation

Install the [Julia](https://julialang.org/) programming language, and do:

```julia-repl
julia> import Pkg; Pkg.add("PDBTools")
```

## Quick Start

```@example index
using PDBTools
atoms = read_pdb(PDBTools.test_dir*"/small.pdb")
```

## What Makes it Special?

- **Simple but Powerful**: Read and write PDB and mmCIF structure files with minimal overhead. 
- **Flexible Atom Selection**: Use intuitive syntax or custom Julia functions.
- **Perfect for MD**: Designed with molecular dynamics workflows in mind.
- **Lightweight**: Focus on atomic data without the overhead of metadata parsing and compact data structures. Can handle very large structures.
- **Performance**: Expect analysis functions to be *fast*. For instance, SASA and hydrogen bonds analyses are among the fastest available.

## Key Features

### Clean Data Structure
Every atom is represented by a simple, accessible structure:

```@example index
atoms[1]
```

```@example index
(name(atoms[1]), resname(atoms[1]), chain(atoms[1]))
```

The use of `InlineStrings` makes the data structure compact, allowing handling millions of atoms is standard computers.

### Intuitive Selection Syntax

Select atoms using simple, readable syntax, similar to that of VMD:
```@example index
selection = select(atoms, "resname ALA and name N CA")
```

### Power of Julia Functions

Leverage Julia's expressiveness for complex selections:
```@example index
selection = select(atoms, atom -> 
    (atom.resname == "ARG" && atom.x < 10) || atom.name == "N"
)
```

Use of Julia functions for selection can also improve performance if dynamic selections are used in critical code.

## See also

`PDBTools.jl` is integrated with [MolSimToolkit.jl](https://github.com/m3g/MolSimToolkit.jl) and [ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl), providing novel and practical tools for molecular dynamics simulations analysis. 

!!! note
    PDBTools prioritizes flexibility over strict format adherence. It's designed for:

    - Molecular dynamics workflows
    - Quick structure analysis
    - Basic PDB/mmCIF file manipulation

    For comprehensive PDB/mmCIF format support, check out [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl) from [BioJulia](https://github.com/BioJulia).
  