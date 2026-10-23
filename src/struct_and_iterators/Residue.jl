"""
    zeta(r::Residue)

Computes the Cα chirality (zeta "virtual" torsion angle - Cα-N-C-Cβ).

Returns the torsion angle or NaN if the residue is not recognized a protein residue or if its a Gly residue.
Expected values are 33.9 ± 3.5 degrees (for one standard deviation). Also see the `zeta_check` function.

# Example

```julia-repl
julia> using PDBTools

julia> protein = select(read_pdb(PDBTools.test_file("structure.pdb")), "protein");

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
        zeta = dihedral(position(r[iCA]), position(r[iN]), position(r[iC]), position(r[iCB]))
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

julia> protein = select(read_pdb(PDBTools.test_file("structure.pdb")), "protein");

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
    protein = select(read_pdb(PDBTools.test_file("structure.pdb")), "protein");
    residues = collect(eachresidue(protein));
    @test zeta(residues[1]) ≈ 33.672016f0
    @test zeta(residues[4]) ≈ 38.015102f0
    @test isnan(zeta(residues[7]))
    @test zeta_check(residues[1]) 
    @test !zeta_check(residues[4]; nsigma=1)
    @test zeta_check(residues[7]) 
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
        resnames = PDBToolsBase.oneletter.(resnames)
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
    ats = read_pdb(PDBTools.test_file("small.pdb"), "protein")
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