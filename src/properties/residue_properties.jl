
const residue_classes = Dict(
    :protein_residues => protein_residues,
    :nucleoside_residues => nucleoside_residues,
)

function _input_residue_classes(residue_classes)
    residue_classes = if isnothing(residue_classes)
        keys(PDBTools.residue_classes)
    elseif residue_classes isa Symbol
        (residue_classes,)
    else
        residue_classes
    end
    if any(!(c in keys(PDBTools.residue_classes)) for c in residue_classes)
        throw(ArgumentError("""\n
            Available residue classes are $(join(':'.*string.(keys(PDBTools.residue_classes)), ", ")) 
            Got: $residue_classes

        """))
    end
    return residue_classes
end

"""
    resname(residue::Union{AbstractString,Char}; residue_classes=nothing)

Returns the residue name, given the one-letter code or residue name. Differently from
`threeletter`, this function will return the force-field name if available in the list
of protein residues.

!!! note
    When retrieving ambiguous codes, protein residue codes will be preferred. Use `class=:nucleoside_residues` to
    retrieve only nucleoside data, or data of any other specific class.

# Examples
```julia-repl
julia> resname("ALA")
"ALA"

julia> resname("GLUP")
"GLUP"

```

"""
function resname(residue::Union{AbstractString,Char}; residue_classes=nothing)
    residue_classes = _input_residue_classes(residue_classes)
    rname, _ = _resname(residue; residue_classes)
    return rname
end

@testitem "resname" begin
    @test resname("ALA") == "ALA"
    @test resname("A") == "ALA"
    @test resname('A') == "ALA"
    @test resname("Alanine") == "ALA"
    @test resname("GLUP") == "GLUP"
    @test resname("HSP") == "HSP"
    @test isnothing(resname("XXX"))
    @test resname("Uridine") == "URD"
    @test resname("A"; residue_classes=:nucleoside_residues) == "ADO"
end

# This variation returns the class where the code was found, for internal processing
function _resname(residue; residue_classes=keys(PDBTools.residue_classes))
    residue_classes = _input_residue_classes(residue_classes)
    code = string(residue)
    # Input: one letter code
    rname = nothing
    cl = nothing
    if length(code) == 1
        rname = nothing
        for class in residue_classes
            rname = findfirst(r -> r.one_letter_code == code, PDBTools.residue_classes[class])
            if !isnothing(rname)
                cl = class
                break
            end
        end
    end
    !isnothing(rname) && return rname, cl
    # Input: some code that matches exactly a three-letter code
    rname = nothing
    cl = nothing
    for class in residue_classes
        _haskey = haskey(PDBTools.residue_classes[class], code)
        if _haskey
            rname = code
            cl = class
            break
        end
    end
    !isnothing(rname) && return rname, cl
    # Input: maybe a full residue name
    rname = nothing
    cl = nothing
    for class in residue_classes
        rname = findfirst(r -> r.name == code, PDBTools.residue_classes[class])
        if !isnothing(rname) 
            cl = class
            break
        end
    end
    !isnothing(rname) && return rname, cl
    return nothing, nothing
end

"""
    threeletter(residue::String; residue_classes=nothing) 

Function to return the three-letter natural-amino acid or nucleoside residue code from the one-letter 
code or residue name. The function is case-insensitive.

!!! note
    When retrieving ambiguous codes, protein residue codes will be preferred. Use `class=:nucleoside_residues` to
    retrieve only nucleoside data, or data of any other specific class.

# Examples

```jldoctest
julia> using PDBTools

julia> threeletter("A")
"ALA"

julia> threeletter("Aspartic acid")
"ASP"

julia> threeletter("HSD")
"HIS"

julia> threeletter("G"; residue_classes=:nucleoside_residues)
"GUO"

```

"""
function threeletter(residue; residue_classes=nothing)
    rname, cl = _resname(residue; residue_classes)
    rname = isnothing(rname) ? nothing : PDBTools.residue_classes[cl][rname].three_letter_code
    return rname
end

@testitem "threeletter" begin
    @test threeletter("HIS") == "HIS"
    @test threeletter("H") == "HIS"
    @test threeletter("E") == "GLU"
    @test threeletter("Histidine") == "HIS"
    @test threeletter("ALA") == "ALA"
    @test threeletter("A") == "ALA"
    @test threeletter('A') == "ALA"
    @test threeletter("Alanine") == "ALA"
    @test threeletter("GLUP") == "GLU"
    @test threeletter("HSP") == "HIS"
    @test isnothing(threeletter("XXX"))
    @test threeletter("G"; residue_classes=:nucleoside_residues) == "GUO"
    @test threeletter("G"; residue_classes=:protein_residues) == "GLY"
    @test_throws "Available residue classes" threeletter("A"; residue_classes=:a)
end

"""
    oneletter(residue::Union{AbstractString,Char})

Function to return a one-letter residue code from the three letter code or residue name. The function is case-insensitive.

### Examples

```julia-repl
julia> oneletter("ALA")
"A"

julia> oneletter("Glutamic acid")
"E"

```

"""
function oneletter(residue::Union{AbstractString,Char})
    code = string(residue)
    rname, cl = _resname(code)
    olc = if isnothing(rname)
        "X"
    else
        residue_classes[cl][rname].one_letter_code
    end
    return olc
end

@testitem "oneletter" begin
    @test oneletter("ALA") == "A"
    @test oneletter("A") == "A"
    @test oneletter('A') == "A"
    @test oneletter("Alanine") == "A"
    @test oneletter("GLUP") == "E"
    @test oneletter("HSP") == "H"
    @test oneletter("XXX") == "X"
    @test oneletter("GUO") == "G"
    @test oneletter("Guanosine") == "G"
end

"""
    residuename(residue::Union{AbstractString,Char})

Returns the long residue name from other residue codes. The function is case-insensitive.

!!! note
    When retrieving ambiguous codes, protein residue codes will be preferred. Use `class=:nucleoside_residues` to
    retrieve only nucleoside data, or data of any other specific class.

### Examples

```julia-repl
julia> residuename("A")
"Alanine"

julia> residuename("Glu")
"Glutamic Acid"

julia> residuename("A"; residue_classes=:nucleoside_residues)
"Adenosine"

""

```

"""
function residuename(residue::Union{AbstractString,Char}; residue_classes=nothing)
    rname, cl = _resname(uppercase(residue); residue_classes)
    residue_name = isnothing(rname) ? nothing : PDBTools.residue_classes[cl][rname].name
    return residue_name
end

@testitem "residuename" begin
    @test residuename("A") == "Alanine"
    @test residuename("A"; residue_classes=:nucleoside_residues) == "Adenosine"
    @test residuename("GLUP") == "Glutamic acid (protonated)"
    @test residuename("Ala") == "Alanine"
end

"""
    Sequence

Wrapper for strings, or vectors of chars, strings, or residue names, to dispatch on 
functions that operate on amino acid sequences.

# Example

```julia-repl
julia> seq = ["Alanine", "Glutamic acid", "Glycine"];

julia> mass(Sequence(seq))
257.2432

julia> seq = "AEG";

julia> mass(Sequence(seq))
257.2432
```

"""
struct Sequence{T}
    s::T
end

"""
    mass(s::Sequence; residue_classes=nothing)

Returns the mass of a sequence of amino acids, given a `Sequence` struct type.

!!! note
    When retrieving ambiguous codes, protein residue codes will be preferred. Use `class=:nucleoside_residues` to
    retrieve only nucleoside data, or data of any other specific class.

# Examples

```julia-repl
julia> seq = ["Alanine", "Glutamic acid", "Glycine"];

julia> mass(Sequence(seq))
257.2432

julia> seq = "AEG";

julia> mass(Sequence(seq))
257.2432

julia> seq = ["ALA", "GLU", "GLY"];

julia> mass(Sequence(seq))
257.2432

julia> mass(Sequence("AUT"))
755.2310170000001
```

"""
function mass(s::Sequence; residue_classes=nothing)
    m = 0.0
    for aa in s.s
        rname, cl = _resname(aa; residue_classes)
        m += PDBTools.residue_classes[cl][rname].mass
    end
    return m
end

@testitem "sequence mass" begin
    @test mass(Sequence("AEG")) ≈ 257.2432
    @test mass(Sequence(["ALA", "GLU", "GLY"])) ≈ 257.2432
    @test mass(Sequence(["A", "E", "G"])) ≈ 257.2432
    @test mass(Sequence(['A', 'E', 'G'])) ≈ 257.2432
    @test mass(Sequence(["Alanine", "Glutamic acid", "Glycine"])) ≈ 257.2432
    @test mass(Sequence("AUT"); residue_classes=:nucleoside_residues) ≈ 755.231017
    @test mass(Sequence(["Adenosine", "Thymidine", "Uridine"])) ≈ 755.231017
end