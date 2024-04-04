#
# List of elements with properties
#
struct Element
    symbol::Symbol
    symbol_string::String
    name::String
    atomic_number::Int
    mass::Float64
    custom::Bool
end

import Base.Broadcast.broadcastable
broadcastable(element::Element) = Ref(element)

#! format: off
const elements = Dict{String,Element}([
    ["X" , "NotFound" ]              .=> Element(:X  ,"X" , "NotFound", 0,  0.00000, false);
    ["H" , "Hydrogen" ]              .=> Element(:H  ,"H" , "Hydrogen", 1,  1.00797, false);
    ["He", "HEL", "Helium" ]         .=> Element(:He ,"He", "Helium", 2,  4.00260, false);
    ["Li", "LI" , "Lithium" ]        .=> Element(:Li ,"Li", "Lithium", 3,  6.941, false);
    ["Be", "BE" , "Beryllium" ]      .=> Element(:Be ,"Be", "Beryllium", 4,  9.01218, false);
    ["B" , "Boron" ]                 .=> Element(:B  ,"B" , "Boron", 5,  10.81, false);
    ["C" , "Carbon" ]                .=> Element(:C  ,"C" , "Carbon", 6,  12.011, false);
    ["N" , "Nitrogen"]               .=> Element(:N  ,"N" , "Nitrogen", 7,  14.0067, false);
    ["O" , "Oxygen" ]                .=> Element(:O  ,"O" , "Oxygen", 8,  15.9994, false);
    ["F" , "Fluorine" ]              .=> Element(:F  ,"F" , "Fluorine", 9,  18.998403, false);
    ["Ne", "NEO", "Neon" ]           .=> Element(:Ne ,"Ne", "Neon", 10,  20.179, false);
    ["Na", "SOD", "NA" , "Sodium" ]  .=> Element(:Na ,"Na", "Sodium", 11,  22.98977, false);
    ["Mg", "MG" , "Magnesium" ]      .=> Element(:Mg ,"Mg", "Magnesium", 12,  24.305, false);
    ["Al", "AL" , "Aluminum" ]       .=> Element(:Al ,"Al", "Aluminum", 13,  26.98154, false);
    ["Si", "SI" , "Silicon" ]        .=> Element(:Si ,"Si", "Silicon", 14,  28.0855, false);
    ["P" , "Phosphorus" ]            .=> Element(:P  ,"P" , "Phosphorus", 15,  30.97376, false);
    ["S" , "Sulfur" ]                .=> Element(:S  ,"S" , "Sulfur", 16,  32.06, false);
    ["Cl", "CL" , "CLA", "Chlorine" ].=> Element(:Cl ,"Cl", "Chlorine", 17,  35.453, false);
    ["Ar", "AR" , "Argon" ]          .=> Element(:Ar ,"Ar", "Argon", 18,  39.948, false);
    ["K" , "POT", "Potassium" ]      .=> Element(:K  ,"K" , "Potassium", 19,  39.0983, false);
    ["Ca", "CAL", "Calcium" ]        .=> Element(:Ca ,"Ca", "Calcium", 20,  40.08, false);
    ["Sc", "SC" , "Scandium" ]       .=> Element(:Sc ,"Sc", "Scandium", 21,  44.9559, false);
    ["Ti", "TI" , "Titanium" ]       .=> Element(:Ti ,"Ti", "Titanium", 22,  47.90, false);
    ["V" , "Vanadium" ]              .=> Element(:V  ,"V" , "Vanadium", 23,  50.9415, false);
    ["Cr", "CR" , "Chromium" ]       .=> Element(:Cr ,"Cr", "Chromium", 24,  51.996, false);
    ["Mn", "MN" , "Manganese" ]      .=> Element(:Mn ,"Mn", "Manganese", 25,  54.9380, false);
    ["Fe", "FE" , "Iron" ]           .=> Element(:Fe ,"Fe", "Iron", 26,  55.847, false);
    ["Co", "CO" , "Cobalt" ]         .=> Element(:Co ,"Co", "Cobalt", 27,  58.9332, false);
    ["Ni", "NI" , "Nickel" ]         .=> Element(:Ni ,"Ni", "Nickel", 28,  58.70, false);
    ["Cu", "CU" , "Copper" ]         .=> Element(:Cu ,"Cu", "Copper", 29,  63.546, false);
    ["Zn", "ZN" , "Zinc" ]           .=> Element(:Zn ,"Zn", "Zinc", 30,  65.38, false);
    ["Ga", "GA" , "Gallium" ]        .=> Element(:Ga ,"Ga", "Gallium", 31,  69.72, false);
    ["Ge", "GE" , "Germanium" ]      .=> Element(:Ge ,"Ge", "Germanium", 32,  72.59, false);
    ["As", "AS" , "Arsenic" ]        .=> Element(:As ,"As", "Arsenic", 33,  74.9216, false);
    ["Se", "SE" , "Selenium" ]       .=> Element(:Se ,"Se", "Selenium", 34,  78.96, false);
    ["Br", "BR" , "Bromine" ]        .=> Element(:Br ,"Br", "Bromine", 35,  79.904, false);
    ["Kr", "KR" , "Krypton" ]        .=> Element(:Kr ,"Kr", "Krypton", 36,  83.80, false);
    ["Rb", "RB" , "Rubidium" ]       .=> Element(:Rb ,"Rb", "Rubidium", 37,  85.4678, false);
    ["Sr", "SR" , "Strontium" ]      .=> Element(:Sr ,"Sr", "Strontium", 38,  87.62, false);
    ["Y" , "Yttrium" ]               .=> Element(:Y  ,"Y" , "Yttrium", 39,  88.9059, false);
    ["Zr", "ZR" , "Zirconium" ]      .=> Element(:Zr ,"Zr", "Zirconium", 40,  91.22, false);
    ["Nb", "NB" , "Niobium" ]        .=> Element(:Nb ,"Nb", "Niobium", 41,  92.9064, false);
    ["Mo", "MO" , "Molybdenum" ]     .=> Element(:Mo ,"Mo", "Molybdenum", 42,  95.94, false);
    ["Tc", "TC" , "Technetium" ]     .=> Element(:Tc ,"Tc", "Technetium", 43,  98.0, false);
    ["Ru", "RU" , "Ruthenium" ]      .=> Element(:Ru ,"Ru", "Ruthenium", 44,  101.07, false);
    ["Rh", "RH" , "Rhodium" ]        .=> Element(:Rh ,"Rh", "Rhodium", 45,  102.9055, false);
    ["Pd", "PD" , "Palladium" ]      .=> Element(:Pd ,"Pd", "Palladium", 46,  106.4, false);
    ["Ag", "AG" , "Silver" ]         .=> Element(:Ag ,"Ag", "Silver", 47,  107.868, false);
    ["Cd", "CAD", "Cadmium" ]        .=> Element(:Cd ,"Cd", "Cadmium", 48,  112.41, false);
    ["In", "IN" , "Indium" ]         .=> Element(:In ,"In", "Indium", 49,  114.82, false);
    ["Sn", "SN" , "Tin" ]            .=> Element(:Sn ,"Sn", "Tin", 50,  118.69, false);
    ["Sb", "SB" , "Antimony" ]       .=> Element(:Sb ,"Sb", "Antimony", 51,  121.75, false);
    ["Te", "TE" , "Tellurium" ]      .=> Element(:Te ,"Te", "Tellurium", 52,  127.60, false);
    ["I" , "Iodine" ]                .=> Element(:I  ,"I" , "Iodine", 53,  126.9045, false);
    ["Xe", "XE" , "Xenon" ]          .=> Element(:Xe ,"Xe", "Xenon", 54,  131.30, false);
    ["Cs", "CES", "Cesium" ]         .=> Element(:Cs ,"Cs", "Cesium", 55,  132.9054, false);
    ["Ba", "BA" , "Barium" ]         .=> Element(:Ba ,"Ba", "Barium", 56,  137.33, false);
    ["La", "LA" , "Lanthanum" ]      .=> Element(:La ,"La", "Lanthanum", 57,  138.9055, false);
    ["Ce", "CER", "Cerium" ]         .=> Element(:Ce ,"Ce", "Cerium", 58,  140.12, false);
    ["Pr", "PR" , "Praseodymium" ]   .=> Element(:Pr ,"Pr", "Praseodymium", 59,  140.9077, false);
    ["Nd", "ND" , "Neodymium" ]      .=> Element(:Nd ,"Nd", "Neodymium", 60,  144.24, false);
    ["Pm", "PM" , "Promethium" ]     .=> Element(:Pm ,"Pm", "Promethium", 61,  145, false);
    ["Sm", "SM" , "Samarium" ]       .=> Element(:Sm ,"Sm", "Samarium", 62,  150.4, false);
    ["Eu", "EU" , "Europium" ]       .=> Element(:Eu ,"Eu", "Europium", 63,  151.96, false);
    ["Gd", "GD" , "Gadolinium" ]     .=> Element(:Gd ,"Gd", "Gadolinium", 64,  157.25, false);
    ["Tb", "TB" , "Terbium" ]        .=> Element(:Tb ,"Tb", "Terbium", 65,  158.9254, false);
    ["Dy", "DY" , "Dysprosium" ]     .=> Element(:Dy ,"Dy", "Dysprosium", 66,  162.50, false);
    ["Ho", "HLM", "Holmium" ]        .=> Element(:Ho ,"Ho", "Holmium", 67,  164.9304, false);
    ["Er", "ER" , "Erbium" ]         .=> Element(:Er ,"Er", "Erbium", 68,  167.26, false);
    ["Tm", "TM" , "Thulium" ]        .=> Element(:Tm ,"Tm", "Thulium", 69,  168.9342, false);
    ["Yb", "YB" , "Ytterbium" ]      .=> Element(:Yb ,"Yb", "Ytterbium", 70,  173.04, false);
    ["Lu", "LU" , "Lutetium" ]       .=> Element(:Lu ,"Lu", "Lutetium", 71,  174.967, false);
    ["Hf", "HAF", "Hafnium" ]        .=> Element(:Hf ,"Hf", "Hafnium", 72,  178.49, false);
    ["Ta", "TA" , "Tantalum" ]       .=> Element(:Ta ,"Ta", "Tantalum", 73,  180.9479, false);
    ["W" , "W"  , "Tungsten" ]       .=> Element(:W  ,"W" , "Tungsten", 74,  183.85, false);
    ["Re", "RE" , "Rhenium" ]        .=> Element(:Re ,"Re", "Rhenium", 75,  186.207, false);
    ["Os", "OS" , "Osmium" ]         .=> Element(:Os ,"Os", "Osmium", 76,  190.2, false);
    ["Ir", "IR" , "Iridium" ]        .=> Element(:Ir ,"Ir", "Iridium", 77,  192.22, false);
    ["Pt", "PT" , "Platinum" ]       .=> Element(:Pt ,"Pt", "Platinum", 78,  195.09, false);
    ["Au", "AU" , "Gold" ]           .=> Element(:Au ,"Au", "Gold", 79,  196.9665, false);
    ["Hg", "MER", "Mercury" ]        .=> Element(:Hg ,"Hg", "Mercury", 80,  200.59, false);
    ["Tl", "TL" , "Thallium" ]       .=> Element(:Tl ,"Tl", "Thallium", 81,  204.37, false);
    ["Pb", "PB" , "Lead" ]           .=> Element(:Pb ,"Pb", "Lead", 82,  207.2, false);
    ["Bi", "BI" , "Bismuth" ]        .=> Element(:Bi ,"Bi", "Bismuth", 83,  208.9804, false);
    ["Po", "PO" , "Polonium" ]       .=> Element(:Po ,"Po", "Polonium", 84,  209, false);
    ["At", "AT" , "Astatine" ]       .=> Element(:At ,"At", "Astatine", 85,  210, false);
    ["Rn", "RN" , "Radon" ]          .=> Element(:Rn ,"Rn", "Radon", 86,  222, false);
    ["Fr", "FR" , "Francium" ]       .=> Element(:Fr ,"Fr", "Francium", 87,  223, false);
    ["Ra", "RA" , "Radium" ]         .=> Element(:Ra ,"Ra", "Radium", 88,  226.0254, false);
    ["Ac", "AC" , "Actinium" ]       .=> Element(:Ac ,"Ac", "Actinium", 89,  227.0278, false);
    ["Th", "TH" , "Thorium" ]        .=> Element(:Th ,"Th", "Thorium", 90,  232.0381, false);
    ["Pa", "PA" , "Protactinium" ]   .=> Element(:Pa ,"Pa", "Protactinium", 91,  231.0359, false);
    ["U" , "U"  , "Uranium" ]        .=> Element(:U  ,"U" , "Uranium", 92,  238.029, false);
])
#! format: on

# Sort element_names by name for faster searching
const element_names = sort(collect(keys(elements)))

""" 
    add_element!(symbol::String, reference_element::PDBTools.Element; elements=PDBTools.elements)

Add a new element to the elements dictionary. If the element already exists, overwrite it.

To remove all custom elements, use `remove_custom_elements!()`.

# Example

```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> using PDBTools

julia> remove_custom_elements!(); # if any

julia> atoms = [ Atom(name="A1"), Atom(name="A2") ];

julia> add_element!("A1", PDBTools.elements["C"])
PDBTools.Element(:C, "C", "Carbon", 6, 12.011, true)

julia> add_element!("A2", PDBTools.elements["N"])
PDBTools.Element(:N, "N", "Nitrogen", 7, 14.0067, true)

julia> element(atoms[1])
"C"

julia> element(atoms[2])
"N"

julia> mass(atoms)
26.017699999999998

julia> remove_custom_elements!(); 
```

Here we repeteadly call `remove_custom_elements!()` to guarantee the proper execution of the
test codes, without any custom elements predefined.

"""
function add_element!(symbol::String, reference_element::Element; elements=elements)
    if symbol in keys(elements)
        @warn """\n
            Element $symbol already exists. Overwriting.

        """ _file=nothing _line=nothing
    end
    elements[symbol] = Element(
            reference_element.symbol,
            reference_element.symbol_string,
            reference_element.name,
            reference_element.atomic_number,
            reference_element.mass,
            true,
        )
    push!(element_names, symbol)
    sort!(element_names)
    return elements[symbol]
end

"""
    remove_custom_elements!()

Remove all custom elements from the elements dictionary.

# Example

```jldoctest
julia> using PDBTools

julia> remove_custom_elements!();

julia> add_element!("GN", PDBTools.elements["N"])
PDBTools.Element(:N, "N", "Nitrogen", 7, 14.0067, true)

julia> element(Atom(name="GN"))
"N"

julia> remove_custom_elements!();

julia> element(Atom(name="GN")) # returns `nothing`

```

Here we repeteadly call `remove_custom_elements!()` to guarantee the proper execution of the
test codes, without any custom elements predefined.

"""
function remove_custom_elements!(elements=PDBTools.elements)
    filter!(name -> !elements[name].custom, element_names)
    filter!(el -> !last(el).custom, PDBTools.elements)
    return PDBTools.elements
end

@testitem "elements" begin
    atoms = readPDB(PDBTools.TESTPDB, "protein")
    @test element.(select(atoms, "residue = 1")) == ["N", "H", "H", "H", "C", "H", "C", "H", "H", "H", "C", "O"]
    @test element_name.(select(atoms, "residue = 1")) == [
        "Nitrogen",
        "Hydrogen",
        "Hydrogen",
        "Hydrogen",
        "Carbon",
        "Hydrogen",
        "Carbon",
        "Hydrogen",
        "Hydrogen",
        "Hydrogen",
        "Carbon",
        "Oxygen",
    ]
    @test atomic_number.(select(atoms, "residue = 1")) == [7, 1, 1, 1, 6, 1, 6, 1, 1, 1, 6, 8]

    # Custom elements
    remove_custom_elements!()
    e_ref = deepcopy(PDBTools.elements)
    add_element!("GN", PDBTools.elements["N"])
    for prop in fieldnames(PDBTools.Element)
        prop == :custom && continue
        @test getfield(PDBTools.elements["GN"],prop) == getfield(PDBTools.elements["N"],prop)
    end
    @test mass(Atom(name="GN")) == mass(Atom(name="N"))
    @test length(PDBTools.element_names) == 268
    remove_custom_elements!()
    @test e_ref == PDBTools.elements
    @test length(PDBTools.element_names) == 267
    @test mass(Atom(name="GN")) === nothing
    remove_custom_elements!()
end
