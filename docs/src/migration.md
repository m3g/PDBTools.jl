# Migrating from PDBTools 3.x to 4.0

PDBTools 4.0 reorganizes the package: the core data structures and the input/output functions
are now defined in [PDBToolsBase.jl](https://github.com/m3g/PDBToolsBase.jl), and re-exported
by PDBTools. For most uses, the functions of the package work as before. This page lists the 
changes that may require updating existing code.

## Periodic boundary conditions: the `unitcell` keyword

This is the change most likely to affect existing code. The `unitcell` keyword argument of 
[`contact_map`](@ref), [`hydrogen_bonds`](@ref), [`sasa_particles`](@ref), [`visualize`](@ref),
and of the transfer free energy functions ([`mvalue`](@ref), `transfer_free_energy`, ...) 
accepts:

| Value                | Meaning                                                      |
|:---------------------|:-------------------------------------------------------------|
| `false` (default)    | No periodic boundary conditions.                             |
| `true`               | Use the unit cell of the [`Structure`](@ref) (read from the file). |
| `UnitCell(...)`      | Use the given unit cell (see [`UnitCell`](@ref)).            |

Matrices, vectors of box sides, and `nothing` are no longer accepted, and throw a `TypeError`.

In PDBTools 3.x:

```julia
atoms = read_pdb("system.pdb")
uc = read_unitcell("system.pdb")
sasa_particles(atoms; unitcell=uc)
contact_map(atoms; unitcell=[85.0, 85.0, 85.0])
hydrogen_bonds(atoms, "protein"; unitcell=nothing)
```

In PDBTools 4.0, the unit cell is read together with the atoms, and can be used directly:

```julia
atoms = read_pdb("system.pdb")
sasa_particles(atoms; unitcell=true)
```

or given explicitly as a `UnitCell`:

```julia
contact_map(atoms; unitcell=UnitCell(85.0, 85.0, 85.0))       # orthorhombic box
contact_map(atoms; unitcell=UnitCell([85.0, 85.0, 85.0]))     # same, from a vector of sides
contact_map(atoms; unitcell=UnitCell(read_unitcell("other.pdb"))) # from a 3×3 matrix
hydrogen_bonds(atoms, "protein")                              # no periodic boundary conditions
```

A `UnitCell` can also be defined from the lattice parameters: `UnitCell(a, b, c, α, β, γ)`, or 
`UnitCell(; a, b, c, α=90, β=90, γ=90)`. It behaves as the 3×3 unit cell matrix, with the 
lattice vectors as columns. [`read_unitcell`](@ref) still returns the matrix.

`unitcell=true` uses the unit cell of the structure given to the function. Since selections 
preserve the unit cell, this also works for subsets of a structure:

```julia
atoms = read_pdb("system.pdb")
protein = select(atoms, "protein")
sasa_particles(protein; unitcell=true) # uses the unit cell of system.pdb
```

If the structure has no unit cell (the file does not define one, or the vector of atoms was
not obtained from a `Structure`), `unitcell=true` throws an `ArgumentError`. The unit cell of a
structure can be set with `atoms.unitcell = UnitCell(...)`.

!!! note
    Files that define the 1Å unit cube in the `CRYST1` (or `_cell`) fields, which indicates 
    structures that are not crystallographic (e.g. NMR or cryo-EM structures), are read 
    without a unit cell.

## Reading functions return a `Structure`

[`read_pdb`](@ref), [`read_mmcif`](@ref) and [`wget`](@ref) return a [`Structure`](@ref), 
instead of a `Vector{Atom}`. A `Structure` is an `AbstractVector{<:Atom}`, and works with all 
functions that accept vectors of atoms. Selections, indexing with ranges, `filter` and `copy` 
of structures return new structures, which preserve the unit cell.

Code that is affected:

- Methods defined for `Vector{<:Atom}` (or `Vector{Atom{Nothing}}`) do not accept structures, 
  and throw a `MethodError`. Define them for `AbstractVector{<:Atom}`:

  ```julia
  f(atoms::Vector{<:Atom}) = ...          # PDBTools 3.x
  f(atoms::AbstractVector{<:Atom}) = ...  # PDBTools 4.0
  ```

- Checks such as `atoms isa Vector` are `false` for structures.

- If a `Vector` is needed, use `collect(atoms)` or `Vector(atoms)` (copies of the vector of 
  atoms), or `atoms.atoms` (the underlying vector, without copying).

- The printed output of structures starts with `Structure{Atom{Nothing}} with N atoms` (and 
  `- with unit cell`, if the unit cell is defined), instead of `Vector{Atom{Nothing}} with N atoms
  with fields:`.

- Files written with [`write_pdb`](@ref) and [`write_mmcif`](@ref) include the unit cell of 
  structures (`CRYST1` or `_cell` fields).

## The `Structure` type

The `Structure` type was experimental in PDBTools 3.x, and was redesigned:

- The metadata is stored in typed fields (currently, only `unitcell`). Arbitrary properties 
  (e.g. `s.filename = "file.pdb"`) are no longer supported.
- The unit cell is a [`UnitCell`](@ref), or `nothing`. Setting `s.unitcell` to a matrix or a 
  vector of sides converts it to a `UnitCell`.
- The type parameters changed to `Structure{T,V}`, where `T` is the type of the atoms and `V` 
  the type of the underlying vector (thus, `Structure{Atom{Nothing}}` matches structures of 
  atoms without custom fields).

## Removed legacy names

The aliases kept for compatibility with older versions were removed:

| Removed                            | Use instead                         |
|:-----------------------------------|:------------------------------------|
| `readPDB(file, ...)`               | [`read_pdb`](@ref)`(file, ...)`     |
| `writePDB(atoms, file, ...)`       | [`write_pdb`](@ref)`(file, atoms, ...)` (note the order of the arguments) |
| `coor(atoms)`                      | [`positions`](@ref)`(atoms)`        |
| `coor(atom)`                       | [`position`](@ref)`(atom)`          |
| `atomic_sasa(atoms, ...)`          | [`sasa_particles`](@ref)`(atoms, ...)` |

## Other changes

- The constants with the paths of test files (`PDBTools.TESTPDB`, `PDBTools.SMALLPDB`, 
  `PDBTools.TESTPBC`, `PDBTools.test_dir`, ...) were removed. The test files are now downloaded
  on demand, with [`PDBTools.test_file`](@ref), which is available after `using LazyArtifacts`: 
  for example, `PDBTools.TESTPDB` is now `PDBTools.test_file("structure.pdb")`, and 
  `PDBTools.SMALLPDB` is `PDBTools.test_file("small.pdb")`.
- PDBTools requires Julia 1.11 or greater.
- The printed output of `SASA` objects (returned by [`sasa_particles`](@ref)) shows only the type of the atomic radii, as in 
  `SASA{StandardAtomicRadii}`.
- [`lattice_to_matrix`](@ref) (and thus [`read_unitcell`](@ref)) returns exact zeros for 
  right angles, which may change the last digits of results obtained with periodic boundary
  conditions.

## For package developers

- The core functions are defined in PDBToolsBase.jl. Packages that only need to read, write,
  and manipulate structures can depend on PDBToolsBase.jl alone, which has lightweight 
  dependencies and is compatible with `juliac --trim`.
- PDBToolsBase.jl replaces MolSimToolkitShared.jl: methods for the shared function names 
  (`distance`, `positions`, `center_of_mass`, `save`, `load`, `coordination_number`, ...) are 
  defined by importing them from PDBToolsBase:

  ```julia
  import MolSimToolkitShared: center_of_mass, wrap   # PDBTools 3.x
  import PDBToolsBase: center_of_mass, wrap          # PDBTools 4.0
  ```
