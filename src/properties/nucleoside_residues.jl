#
# Data for natural nucleosides
#
@kwdef struct NucleosideResidue
    name::String
    three_letter_code::String
    one_letter_code::String
    type::String
    mono_isotopic_mass::Float64
    mass::Float64
    charge::Int
    custom::Bool = false
end

#! format: off
const nucleoside_residues = OrderedDict{String,NucleosideResidue}(
    "ADO" => NucleosideResidue("Adenosine", "ADO", "A", "Purine",     267.094691, 267.094691, 0, false),
    "GUO" => NucleosideResidue("Guanosine", "GUO", "G", "Purine",     283.090169, 283.090169, 0, false),
    "URD" => NucleosideResidue("Uridine",   "URD", "U", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "CYD" => NucleosideResidue("Cytidine",  "CYD", "C", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "THD" => NucleosideResidue("Thymidine", "THD", "T", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "INO" => NucleosideResidue("Inosine",   "INO", "I", "Purine",     268.090169, 268.090169, 0, false),
)
#! format: on

"""
    add_nucleoside_residue!(resname::String, reference_residue::PDBTools.NucleosideResidue)

Add a custom nucleoside residue to the list of nucleoside residues. Returns the added
`NucleosideResidue`. Use `remove_custom_nucleoside_residues!()` to undo.

# Example

```jldoctest
julia> using PDBTools

julia> remove_custom_nucleoside_residues!();

julia> add_nucleoside_residue!("MYN", PDBTools.nucleoside_residues["ADO"])
PDBTools.NucleosideResidue("MYN", "ADO", "A", "Purine", 267.094691, 267.094691, 0, true)

julia> isnucleoside(Atom(resname="MYN"))
true

julia> remove_custom_nucleoside_residues!(); # clean up
```

"""
function add_nucleoside_residue!(resname::String, reference_residue::PDBTools.NucleosideResidue)
    if haskey(PDBTools.nucleoside_residues, resname)
        @warn """\n
            Residue $resname already exists in the list of nucleoside residues. Overwriting.

        """ _file = nothing _line = nothing
    end
    PDBTools.nucleoside_residues[resname] = PDBTools.NucleosideResidue(
        resname,
        reference_residue.three_letter_code,
        reference_residue.one_letter_code,
        reference_residue.type,
        reference_residue.mono_isotopic_mass,
        reference_residue.mass,
        reference_residue.charge,
        true,
    )
    return PDBTools.nucleoside_residues[resname]
end

"""
    remove_custom_nucleoside_residues!()

Remove all custom nucleoside residues from the nucleoside residue list.

# Example

```jldoctest
julia> using PDBTools

julia> remove_custom_nucleoside_residues!();

julia> add_nucleoside_residue!("MYN", PDBTools.nucleoside_residues["ADO"])
PDBTools.NucleosideResidue("MYN", "ADO", "A", "Purine", 267.094691, 267.094691, 0, true)

julia> isnucleoside(Atom(resname="MYN"))
true

julia> remove_custom_nucleoside_residues!();

julia> isnucleoside(Atom(resname="MYN"))
false
```

"""
remove_custom_nucleoside_residues!() = filter!(r -> !last(r).custom, nucleoside_residues)

@testitem "nucleoside residues" begin
    using PDBTools
    remove_custom_nucleoside_residues!()

    @test isnucleoside(Atom(resname="ADO"))
    @test isnucleoside(Atom(resname="GUO"))
    @test isnucleoside(Atom(resname="URD"))
    @test isnucleoside(Atom(resname="CYD"))
    @test isnucleoside(Atom(resname="THD"))
    @test isnucleoside(Atom(resname="INO"))
    @test !isnucleoside(Atom(resname="ALA"))

    @test ispurine(Atom(resname="ADO"))
    @test ispurine(Atom(resname="GUO"))
    @test ispurine(Atom(resname="INO"))
    @test !ispurine(Atom(resname="URD"))

    @test ispyrimidine(Atom(resname="URD"))
    @test ispyrimidine(Atom(resname="CYD"))
    @test ispyrimidine(Atom(resname="THD"))
    @test !ispyrimidine(Atom(resname="ADO"))

    # custom residue add/remove
    add_nucleoside_residue!("MYN", PDBTools.nucleoside_residues["ADO"])
    @test isnucleoside(Atom(resname="MYN"))
    @test ispurine(Atom(resname="MYN"))
    remove_custom_nucleoside_residues!()
    @test !isnucleoside(Atom(resname="MYN"))

    # string selection
    pdb = read_pdb(PDBTools.TESTPDB)
    @test length(select(pdb, "nucleoside")) == 0  # no nucleosides in protein test structure
end