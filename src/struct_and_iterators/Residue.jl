"""
    Residue

Residue data structure. 

The Residue structure carries the properties of the residue or molecule of the atoms
it contains, but it does not copy the original vector of atoms, only the residue
meta data for each residue. Thus, changes in the residue atoms will be reflected in the
original vector of atoms.

### Example

```jldoctest
julia> using PDBTools

julia> pdb = wget("1LBD");

julia> residues = collect(eachresidue(pdb))
238-element Vector{Residue}[
    SER225A
    ALA226A
    ⋮
    MET461A
    THR462A
]

julia> resnum.(residues[1:3])
3-element Vector{Int32}:
 225
 226
 227

julia> residues[5].chain
"A"

julia> residues[8].range
52:58

julia> mass(residues[1])
82.0385f0

```

"""
struct Residue{T<:Atom,Vec<:AbstractVector{T}} <: AbstractStructuralElement{T}
    atoms::Vec
    range::UnitRange{Int}
    resname::StringType
    residue::Int32
    resnum::Int32
    chain::StringType
    model::Int32
    segname::StringType
end

# Necessary for the interface: define the same_struct_element function
same_struct_element(::Type{Residue}, at1::Atom, at2::Atom) = same_residue(at1, at2)
Base.in(at::Atom, r::Residue) = same_struct_element(Residue, at, first(r))

# Constructors
function Residue(atoms::AbstractVector{<:Atom}, range::AbstractRange{<:Integer})
    _check_unique(Residue, atoms, range)
    at1 = first(atoms[range])
    Residue(atoms, UnitRange{Int}(range), resname(at1), residue(at1), resnum(at1), chain(at1), model(at1), segname(at1))
end
Residue(atoms::AbstractVector{<:Atom}) = Residue(atoms, 1:length(atoms))

"""
    zeta(r::Residue)

Computes the Cα chirality (zeta "virtual" torsion angle - Cα-N-C-Cβ).

Returns the torsion angle or NaN if the residue is not recognized a protein residue or if its a Gly residue.
Expected values are 33.9 ± 3.5 degrees (for one standard deviation). Also see the `zeta_check` function.

# Example

```julia-repl
julia> using PDBTools

julia> protein = select(read_pdb(PDBTools.TESTPDB), "protein");

julia> residues = collect(eachresidue(protein));

julia> zeta(residues[1])
33.67202f0
```

"""
function zeta(r::Residue)
    if isprotein(r) && resname(r) != "GLY"
        iCA = findfirst(at -> name(at) == "CA", r)
        iN = findfirst(at -> name(at) == "N", r)
        iC = findfirst(at -> name(at) == "C", r)
        iCB = findfirst(at -> name(at) == "CB", r)
        zeta = dihedral(coor(r[iCA]), coor(r[iN]), coor(r[iC]), coor(r[iCB]))
        return zeta
    else
        return NaN
    end
end

"""
    zeta_check(r::Residue; nsigma=2)

Checks if the Cα chirality falls into expected ranges. See the `zeta` function for further
information. The expected mean is 33.9 degrees with a standard deviation of 3.5 degrees.
By default, `nsigma=2`, implying that the function returns `true` if the torsion falls 
within two standard deviations from the mean.

See: https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/manual/manappa.html

# Example

```jldoctest
julia> using PDBTools

julia> protein = select(read_pdb(PDBTools.TESTPDB), "protein");

julia> residues = collect(eachresidue(protein));

julia> zeta_check(residues[1])
true
```

"""
function zeta_check(r::Residue; nsigma=2)
    # data from: https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/manual/manappa.html
    mean = 33.9 # degrees 
    margin = nsigma * 3.5 # std
    if isprotein(r) && resname(r) != "GLY"
        return mean - margin < zeta(r) < mean + margin
    else
        return true
    end
end

@testitem "zeta" begin
    using PDBTools
    protein = select(read_pdb(PDBTools.TESTPDB), "protein");
    residues = collect(eachresidue(protein));
    @test zeta(residues[1]) ≈ 33.672016f0
    @test zeta(residues[4]) ≈ 38.015102f0
    @test isnan(zeta(residues[7]))
    @test zeta_check(residues[1]) 
    @test !zeta_check(residues[4]; nsigma=1)
    @test zeta_check(residues[7]) 
end


"""
    eachresidue(atoms::AbstractVector{<:Atom})

Iterator for the residues (or molecules) of a selection. 

### Example

```jldoctest
julia> using PDBTools

julia> atoms = wget("1LBD");

julia> eachresidue(atoms)
 Residue iterator with length = 238

julia> collect(eachresidue(atoms))
238-element Vector{Residue}[
    SER225A
    ALA226A
    ⋮
    MET461A
    THR462A
]
```

"""
eachresidue(atoms::AbstractVector{<:Atom}) = EachStructuralElement{Residue}(atoms)

# Specific getters for this type

"""
    name(residue::Residue)

Returns the name (residue name) of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `String`: The residue name.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> name(residues[1])
"ALA"
```

"""
name(residue::Residue) = residue.resname

"""
    resname(residue::Residue)

Returns the residue name of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `String`: The residue name.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> resname(residues[1])
"ALA"
```

"""
resname(residue::Residue) = residue.resname

"""
    residue(residue::Residue)

Returns the sequential residue (molecule) number of the residue in the structure.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `Int32`: The sequential residue number.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> residue(residues[1])
1
```

"""
residue(residue::Residue) = residue.residue

"""
    resnum(residue::Residue)

Returns the residue number as written in the PDB file.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `Int32`: The residue number from the PDB file.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> resnum(residues[1])
1
```

"""
resnum(residue::Residue) = residue.resnum

"""
    chain(residue::Residue)

Returns the chain identifier of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `String`: The chain identifier.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> chain(residues[1])
"A"
```

"""
chain(residue::Residue) = residue.chain

"""
    model(residue::Residue)

Returns the model number of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `Int32`: The model number.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> model(residues[1])
1
```

"""
model(residue::Residue) = residue.model

"""
    segname(residue::Residue)

Returns the segment name of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `String`: The segment name.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> segname(residues[1])
"PROT"
```

"""
segname(residue::Residue) = residue.segname

"""
    mass(residue::Residue)

Returns the total mass of all atoms in the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `Float32`: The total mass in atomic mass units (amu).

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> mass(residues[1])
82.0385f0
```

"""
mass(residue::Residue) = mass(@view residue.atoms[residue.range])

"""
    get_atoms(residue::Residue)

Returns a view of the atoms that belong to the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `SubArray{Atom}`: A view of the atoms in the residue.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> res_atoms = get_atoms(residues[1]);

julia> length(res_atoms)
12
```

"""
get_atoms(residue::Residue) = @view(residue.atoms[residue.range])

"""
    charge(residue::Residue)

Returns the formal charge of the residue. For protein residues, the charge is obtained
from the standard protein residue database. For non-protein residues, a warning is issued
and 0 is returned.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `Int`: The formal charge of the residue.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.TESTPDB, "protein and residue 1");

julia> residues = collect(eachresidue(atoms));

julia> charge(residues[1])
1
```

"""
function charge(residue::Residue)  
    if residue.resname in keys(protein_residues)
        return protein_residues[residue.resname].charge
    else
        @warn "Residue $(residue.resname) not found in protein residues. Assuming residue charge = 0"
        return 0
    end
end

@testitem "Residue iterator" begin
    using PDBTools
    atoms = read_pdb(PDBTools.TESTPDB, "protein")
    residues = eachresidue(atoms)
    @test length(residues) == 104
    @test name(Residue(atoms, 1:12)) == "ALA"
    @test resname(Residue(atoms, 1:12)) == "ALA"
    @test Residue(atoms, 1:12).range == 1:12
    @test Residue(atoms[1:12]).range == 1:12
    @test firstindex(residues) == 1
    @test lastindex(residues) == 104
    @test resname(last(residues)) == "THR"
    @test residue(last(residues)) == 104
    @test chain(last(residues)) == "A"
    @test_throws ArgumentError residues[1]
    residues = collect(eachresidue(atoms))
    @test atoms[1] in residues[1]
    @test !(atoms[1] in residues[2])
    @test atoms[1458] in residues[104]
    @test index.(filter(at -> name(at) in ("N", "HG1"), residues[2])) == [13, 21]
    @test length(get_atoms(residues[1])) == 12
    @test findall(at -> name(at) in ("N", "HG1"), residues[2]) == [1, 9]
    @test_throws ArgumentError Residue(atoms, 1:15)
end

#
# io show functions
#
function Base.show(io::IO, residue::Residue)
    compact = get(io, :compact, false)::Bool
    if compact
        print(io, "$(name(residue))$(resnum(residue))$(chain(residue))")
    else
        println(io, " Residue of name $(name(residue)) with $(length(residue)) atoms.")
        show(IOContext(io, :type => false), @view residue.atoms[residue.range])
    end
end

@testitem "residue show" begin
    using PDBTools
    using ShowMethodTesting
    ENV["LINES"] = 10
    ENV["COLUMNS"] = 120
    ats = read_pdb(PDBTools.SMALLPDB)
    r = eachresidue(ats)
    @test parse_show(r; repl=Dict("PDBTools." => "")) ≈ "Residue iterator with length = 3"
    rc = collect(r)
    @test parse_show(rc; repl=Dict("PDBTools." => "")) ≈ """
    3-element Vector{Residue}[ 
        ALA1A
        CYS2A
        ASP3A
    ]
    """
    @test parse_show(rc[1]; repl=Dict(r"^((?:[^\n]*\n){3}).*"s => s"\1")) ≈ """
     Residue of name ALA with 12 atoms.
    index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
    """
end

#
# Properties of residues
#
# Deal with alternate conformers with the second test
#
isprotein(residue::Residue) = 
    haskey(protein_residues, resname(residue)) || 
    (length(resname(residue)) == 4 && haskey(protein_residues, @view(resname(residue)[2:4])))

export isprotein
export isacidic, isaliphatic, isaromatic, isbasic, ischarged
export ishydrophobic, isneutral, isnonpolar, ispolar
export iswater

isacidic(r::Residue) = isprotein(r) && protein_residues[r.resname].type == "Acidic"
isaliphatic(r::Residue) = isprotein(r) && protein_residues[r.resname].type == "Aliphatic"
isaromatic(r::Residue) = isprotein(r) && protein_residues[r.resname].type == "Aromatic"
isbasic(r::Residue) = isprotein(r) && protein_residues[r.resname].type == "Basic"
ischarged(r::Residue) = isprotein(r) && protein_residues[r.resname].charge != 0
isneutral(r::Residue) = isprotein(r) && protein_residues[r.resname].charge == 0
ishydrophobic(r::Residue) = isprotein(r) && protein_residues[r.resname].hydrophobic
ispolar(r::Residue) = isprotein(r) && protein_residues[r.resname].polar
isnonpolar(r::Residue) = isprotein(r) && !ispolar(r)

isacidic(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].type == "Acidic"
isaliphatic(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].type == "Aliphatic"
isaromatic(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].type == "Aromatic"
isbasic(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].type == "Basic"
ischarged(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].charge != 0
isneutral(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].charge == 0
ishydrophobic(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].hydrophobic
ispolar(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].polar
isnonpolar(atom::Atom) = isprotein(atom) && !ispolar(atom)

const water_residues = ["HOH", "OH2", "TIP3", "TIP3P", "TIP4P", "TIP5P", "TIP7P", "SPC", "SPCE"]
iswater(r::Residue; water_residues=water_residues) = r.resname in water_residues
iswater(atom::Atom; water_residues=water_residues) = atom.resname in water_residues

@testitem "residue of atom" begin
    pdb = read_pdb(PDBTools.TESTPDB)
    glu = select(pdb, "resname GLU")
    @test isacidic(glu[1])
    @test !isaliphatic(glu[1])
    @test !isaromatic(glu[1])
    @test !isbasic(glu[1])
    @test ischarged(glu[1])
    @test !isneutral(glu[1])
    @test !ishydrophobic(glu[1])
    @test ispolar(glu[1])
    @test !isnonpolar(glu[1])
    @test !iswater(glu[1])
    phe = select(pdb, "resname PHE")
    @test !isacidic(phe[1])
    @test !isaliphatic(phe[1])
    @test isaromatic(phe[1])
    @test !isbasic(phe[1])
    @test !ischarged(phe[1])
    @test isneutral(phe[1])
    @test ishydrophobic(phe[1])
    @test !ispolar(phe[1])
    @test isnonpolar(phe[1])
    @test !iswater(glu[1])
    wat = select(pdb, "water")
    @test iswater(wat[1])
    pdb = [ Atom(resname="AGLY"), Atom(resname="BGLY")]
    @test length(eachresidue(pdb)) == 1
    @test resname.(eachresidue(pdb)) == ["AGLY"] 
end

@testitem "full residue" begin
    pdb = read_pdb(PDBTools.TESTPDB)
    glu_atoms = select(pdb, "resname GLU")
    glu = collect(eachresidue(glu_atoms))
    @test charge(glu[1]) == -1
    @test isacidic(glu[1])
    @test !isaliphatic(glu[1])
    @test !isaromatic(glu[1])
    @test !isbasic(glu[1])
    @test ischarged(glu[1])
    @test !isneutral(glu[1])
    @test !ishydrophobic(glu[1])
    @test ispolar(glu[1])
    @test !isnonpolar(glu[1])
    @test !iswater(glu[1])
    @test mass(glu[1]) ≈ 128.1077 atol = 1e-3
    @test mass(glu[begin]) ≈ mass(glu[end])
    @test segname(glu[begin]) == "PROT"
    @test segname(glu[end]) == "PROT"
    phe_atoms = select(pdb, "resname PHE")
    phe = collect(eachresidue(phe_atoms))
    @test charge(phe[1]) == 0
    @test !isacidic(phe[1])
    @test !isaliphatic(phe[1])
    @test isaromatic(phe[1])
    @test !isbasic(phe[1])
    @test !ischarged(phe[1])
    @test isneutral(phe[1])
    @test ishydrophobic(phe[1])
    @test !ispolar(phe[1])
    @test isnonpolar(phe[1])
    @test !iswater(glu[1])
    water = select(pdb, "water and resnum <= 3")
    watresiter = eachresidue(water)
    watres = collect(eachresidue(water))
    @test iswater(watres[1])
    @test charge(watres[1]) == 0
end

"""
    residue_ticks(
        atoms (or) residues (or) residue iterator; 
        first=nothing, last=nothing, stride=1, oneletter=true, serial=false
    )

Returns a tuple with residue numbers and residue names for the given atoms, to be used as tick labels in plots.

The structure data can be provided a vector of `Atom`s, a vector of `Residue`s or an `eachresidue` iterator. 

`first` and `last` optional keyword parameters are integers that refer to the residue numbers to be included. 
The `stride` option can be used to skip residues and declutter the tick labels.

If `oneletter` is `false`, three-letter residue codes are returned. Residues with unknown names will be 
named `X` or `XXX`. 

If `serial=true` the sequential residue index will be used as the index of the ticks. If instead
`serial=false`, the positions will be set to the residue numbers.

# Examples

```jldoctest
julia> using PDBTools

julia> atoms = wget("1LBD", "protein");

julia> residue_ticks(atoms; stride=50) # Vector{<:Atom} as input
(Int32[225, 275, 325, 375, 425], ["S225", "Q275", "L325", "L375", "L425"])

julia> residue_ticks(atoms; first=235, last=240) # first=10
(Int32[235, 236, 237, 238, 239, 240], ["I235", "L236", "E237", "A238", "E239", "L240"])

julia> residue_ticks(eachresidue(atoms); stride=50) # residue iterator as input
(Int32[225, 275, 325, 375, 425], ["S225", "Q275", "L325", "L375", "L425"])

julia> residue_ticks(collect(eachresidue(atoms)); stride=50) # Vector{Residue} as input
(Int32[225, 275, 325, 375, 425], ["S225", "Q275", "L325", "L375", "L425"])

julia> residue_ticks(atoms; first=10, stride=50, serial=true) # using serial=true
(10:50:210, ["R234", "K284", "R334", "S384", "E434"])

```

The resulting tuple of residue numbers and labels can be used as `xticks` in `Plots.plot`, for example.

"""
function residue_ticks(residues::Union{AbstractVector{Residue},EachStructuralElement{<:Residue}};
    first=nothing, last=nothing, stride=1,
    oneletter::Bool=true,
    serial::Bool=false,
)
    resnames = resname.(residues)
    if oneletter
        resnames = PDBTools.oneletter.(resnames)
    end
    resnums = resnum.(residues)
    ticklabels = resnames .* string.(resnums)
    if !serial && (!allunique(resnums) || !issorted(resnums))
        @warn """\n
            Residue numbers are not unique and/or are not sorted. Using serial indexing (1:$stride:$(lastindex(residues))) to define x-tick positions.

        """ _file = nothing _line = nothing
        serial = true
    end
    if serial
        first = isnothing(first) ? 1 : first
        last = isnothing(last) ? length(residues) : last
        if first < firstindex(resnums) || last > lastindex(resnums)
            throw(ArgumentError("""\n
                First or last residue index out of residue number range: $(firstindex(resnums)) to $(lastindex(resnums)).

            """))
        end
        ticks = (first:stride:last, ticklabels[first:stride:last])
    else
        first = isnothing(first) ? 1 : findfirst(==(first), resnums)
        last = isnothing(last) ? length(residues) : findfirst(==(last), resnums)
        if (isnothing(first) || isnothing(last))
            throw(ArgumentError("""\n
                First or last residue index out of residue number range: $(minimum(resnums)) to $(maximum(resnums)).

            """))
        end
        ticks = (resnums[first:stride:last], ticklabels[first:stride:last])
    end
    return ticks
end
function residue_ticks(atoms::AbstractVector{<:Atom}; kargs...)
    residues = eachresidue(atoms)
    return residue_ticks(residues; kargs...)
end

@testitem "residue_ticks" begin
    # residue indices start with 1
    atoms = wget("1UBQ")
    @test residue_ticks(eachresidue(atoms), stride=20, serial=false) == ([1, 21, 41, 61, 81, 101, 121], ["M1", "D21", "Q41", "I61", "X81", "X101", "X121"])
    @test residue_ticks(collect(eachresidue(atoms)), stride=20, serial=false) == ([1, 21, 41, 61, 81, 101, 121], ["M1", "D21", "Q41", "I61", "X81", "X101", "X121"])
    @test residue_ticks(atoms, stride=20, serial=false) == ([1, 21, 41, 61, 81, 101, 121], ["M1", "D21", "Q41", "I61", "X81", "X101", "X121"])
    atoms = select(atoms, "protein")
    # serial indexing
    @test residue_ticks(atoms; first=10, stride=10, serial=true) == (10:10:70, ["G10", "S20", "I30", "Q40", "L50", "N60", "V70"])
    @test residue_ticks(atoms; first=1, stride=10, serial=true) == (1:10:71, ["M1", "K11", "D21", "Q31", "Q41", "E51", "I61", "L71"])
    @test residue_ticks(atoms; first=10, stride=1, last=15, serial=true) == (10:1:15, ["G10", "K11", "T12", "I13", "T14", "L15"])
    @test residue_ticks(atoms; first=10, stride=2, last=15, serial=true) == (10:2:14, ["G10", "T12", "T14"])
    # resnum indexing
    @test residue_ticks(atoms; stride=20, serial=false) == ([1, 21, 41, 61], ["M1", "D21", "Q41", "I61"])
    @test residue_ticks(atoms; stride=20, first=2, serial=false) == ([2, 22, 42, 62], ["Q2", "T22", "R42", "Q62"])
    @test residue_ticks(atoms; stride=20, last=42, serial=false) == ([1, 21, 41], ["M1", "D21", "Q41"])
    @test residue_ticks(atoms; stride=20, last=42, first=2, serial=false) == ([2, 22, 42], ["Q2", "T22", "R42"])
    # Shift resnums
    for at in atoms
        at.resnum += 10
    end
    @test residue_ticks(atoms; stride=20, serial=false) == ([11, 31, 51, 71], ["M11", "D31", "Q51", "I71"])
    @test_throws ArgumentError residue_ticks(atoms; stride=20, first=2, serial=false)
    @test residue_ticks(atoms; stride=20, first=13, serial=false) == ([13, 33, 53, 73], ["I13", "I33", "L53", "K73"])
    @test residue_ticks(atoms; stride=20, last=42, serial=false) == ([11, 31], ["M11", "D31"])
    @test_throws ArgumentError residue_ticks(atoms; stride=20, first=42, last=90, serial=false)
    @test residue_ticks(atoms; stride=20, first=42, last=85, serial=false) == ([42, 62, 82], ["D42", "D62", "R82"])
    # serial indexing
    @test residue_ticks(atoms; stride=20, serial=true) == (1:20:61, ["M11", "D31", "Q51", "I71"])
    @test_throws ArgumentError residue_ticks(atoms; stride=20, first=0, serial=true)
    @test residue_ticks(atoms; stride=20, first=13, serial=true) == (13:20:73, ["I23", "K43", "G63", "L83"])
    @test residue_ticks(atoms; stride=20, last=42, serial=true) == (1:20:41, ["M11", "D31", "Q51"])
    @test_throws ArgumentError residue_ticks(atoms; stride=20, first=42, last=90, serial=true)
    @test residue_ticks(atoms; stride=20, first=32, last=75, serial=true) == (32:20:72, ["D42", "D62", "R82"])
    # residue indices do not start with 1
    atoms = wget("1LBD", "protein")
    @test residue_ticks(atoms, stride=38, serial=false) == ([225, 263, 301, 339, 377, 415, 453], ["S225", "D263", "L301", "S339", "N377", "F415", "E453"])
    @test residue_ticks(atoms; stride=1, first=227, last=231, serial=false) == ([227, 228, 229, 230, 231], ["N227", "E228", "D229", "M230", "P231"])
    @test residue_ticks(atoms; stride=2, first=227, last=231, serial=false) == ([227, 229, 231], ["N227", "D229", "P231"])
    @test residue_ticks(atoms; stride=2, first=227, last=231, serial=true) == (227:2:231, ["L451", "E453", "L455"])
    # three-letter return codes
    @test residue_ticks(atoms; stride=2, first=227, last=231, oneletter=false, serial=false) == ([227, 229, 231], ["ASN227", "ASP229", "PRO231"])

    # non-unique residue numbers
    ats = read_pdb(PDBTools.SMALLPDB, "protein")
    ats_dup = vcat(ats, copy(ats))
    @test residue_ticks(ats_dup) == (1:1:6, ["A1", "C2", "D3", "A1", "C2", "D3"])
    @test residue_ticks(ats_dup; stride=3) == (1:3:4, ["A1", "A1"])
    @test residue_ticks(ats_dup; first=2, last=5) == (2:1:5, ["C2", "D3", "A1", "C2"])
    @test residue_ticks(ats_dup; serial=false) == (1:1:6, ["A1", "C2", "D3", "A1", "C2", "D3"])

    # additional input errors
    @test_throws ArgumentError residue_ticks(ats; first=0)
    @test_throws ArgumentError residue_ticks(ats; last=36)
    @test_throws ArgumentError residue_ticks(ats; serial=true, first=0)
    @test_throws ArgumentError residue_ticks(ats; serial=true, last=36)

    # resnums out of order
    for residue in eachresidue(ats)
        if resnum(residue) == 2
            for at in residue
                at.resnum = 4
            end
        end
    end
    @test residue_ticks(ats; serial=false) == (1:1:3, ["A1", "C4", "D3"])
    @test residue_ticks(ats; serial=true) == (1:1:3, ["A1", "C4", "D3"])

end