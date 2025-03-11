export add_hydrogens!

"""
    add_hydrogens!(atoms::AbstractVector{<:Atom}; pH=7.0, obabel="obabel", debug=false)

Add hydrogens to a PDB file using Open Babel. 

# Arguments

- `atoms::AbstractVector{<:Atom}`: structure (usually PDB file of a protein) to add hydrogens to.
- `pH`: the pH of the solution. Default is 7.0.
- `obabel`: path to the obabel executable. Default is "obabel".
- `debug`: if true, print the output message from obabel. Default is false.

!!! note
    This function requires the installation of [OpenBabel](http://openbabel.org/).
    Please cite the corresponding reference if using it.

# Example

```julia-repl
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.TESTPDB, "protein and not element H");

julia> add_hydrogens!(atoms)
   Vector{Atom{Nothing}} with 1459 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  1.00  0.00     1       -         1
       2   CA     ALA     A        1        1   -8.483  -14.912   -6.726  1.00  0.00     1       -         2
       3   CB     ALA     A        1        1   -9.383  -14.465   -7.880  1.00  0.00     1       -         3
                                                       ⋮ 
    1457    H     THR     A      104      208    5.886  -10.722   -7.797  1.00  0.00     1       -      1457
    1458    H     THR     A      104      208    5.871  -10.612   -9.541  1.00  0.00     1       -      1458
    1459    H     THR     A      104      208    6.423  -12.076   -8.762  1.00  0.00     1       -      1459
```

"""
function add_hydrogens!(
    atoms::AbstractVector{<:Atom};
    pH::Real=7.0,
    obabel="obabel",
    debug=false,
)
    tmpfile = tempname() * ".pdb"
    tmpfile_out = tempname() * ".pdb"
    write_pdb(tmpfile, atoms)
    if isnothing(Sys.which(obabel))
        throw(ArgumentError("""

            "$obabel" executable not found.
            
            Is it installed? Provide the path to the executable with the `obabel` keyword.

            Run with debug=true to see the error message from obabel.
        """))
    end
    try
        if debug
            readchomp(pipeline(`$obabel $tmpfile -O $tmpfile_out -p $pH`))
        else
            readchomp(pipeline(`$obabel $tmpfile -O $tmpfile_out -p $pH`; stderr=devnull))
        end
    catch
        throw(ArgumentError("""
            Error running obabel. 
            Run with debug=true to see the error message from obabel.
        """))
    end
    atoms_read = read_pdb(tmpfile_out)
    sort!(atoms_read; by=at -> resnum(at))
    setproperty!.(atoms_read, :index, eachindex(atoms_read))
    atoms .= atoms_read[1:length(atoms)]
    append!(atoms, atoms_read[length(atoms)+1:end])
    return atoms
end

@testitem "add_hydrogens!" begin
    using PDBTools
    atoms = read_pdb(PDBTools.TESTPDB, "protein and not element H")
    @test length(atoms) == 781
    if !isnothing(Sys.which("obabel"))
        add_hydrogens!(atoms)
        @test length(atoms) == 1459
        @test count(sel"element H", atoms) == 678
    end
    @test_throws TypeError add_hydrogens!(atoms; pH="A")
    @test_throws ArgumentError add_hydrogens!(atoms; obabel="nonexistent")
end

"""
    center_of_mass(atoms::AbstractVector{<:Atom})

Calculate the center of mass of the atoms.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> center_of_mass(atoms)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  -5.584422772707132
 -13.110413081059928
  -7.139970851058855
```

"""
function center_of_mass(atoms::AbstractVector{<:Atom})
    totmass = mass(atoms)
    center = SVector(
        sum(at.x*mass(at) for at in atoms), 
        sum(at.y*mass(at) for at in atoms), 
        sum(at.z*mass(at) for at in atoms)
    ) ./ totmass
    return center
end

@testitem "center_of_mass" begin
    using PDBTools
    atoms = [ Atom(name="C", x=-1.0, y=-1.0, z=-1.0), Atom(name="C", x=1.0, y=1.0, z=1.0) ]
    @test center_of_mass(atoms) ≈ [0.0, 0.0, 0.0] atol = 1e-10
    atoms = read_pdb(PDBTools.SMALLPDB)
    @test center_of_mass(atoms) ≈ [-5.584422752942997, -13.110413157869903, -7.139970815730879] atol = 1e-7
end

"""
    moveto!(atoms::AbstractVector{<:Atom}; center::AbstractVector{<:Real}=SVector(0.0, 0.0, 0.0))

Move the center of mass of the atoms to the specified `center` position, which defaults to the origin.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> center_of_mass(atoms)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  -5.584422772707132
 -13.110413081059928
  -7.139970851058855

julia> moveto!(atoms; center = [1.0, 2.0, 3.0]);

julia> center_of_mass(atoms)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 1.0000000263619948
 1.9999999934852166
 2.9999999509668918
```
"""
function moveto!(atoms::AbstractVector{<:Atom}; center::AbstractVector{<:Real}=SVector(0.0, 0.0, 0.0))
    cm = center_of_mass(atoms)
    for at in atoms
        at.x += center[1] - cm[1]
        at.y += center[2] - cm[2]
        at.z += center[3] - cm[3]
    end
    return atoms
end

@testitem "moveto!" begin
    using PDBTools
    atoms = read_pdb(PDBTools.SMALLPDB)
    @test center_of_mass(moveto!(atoms)) ≈ [0.0, 0.0, 0.0] atol = 1e-7
    @test center_of_mass(moveto!(atoms; center=[1.0, 2.0, 3.0])) ≈ [1.0, 2.0, 3.0] atol=1e-7
end