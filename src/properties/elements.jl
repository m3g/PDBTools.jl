#
# List of elements with properties
#
# vdw radius from: https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page) 
#
@kwdef struct Element
    symbol::Symbol
    symbol_string::String7
    name::String
    atomic_number::Int
    mass::Float32
    custom::Bool
    vdw_radius::Float32
end

import Base.Broadcast.broadcastable
broadcastable(element::Element) = Ref(element)

#! format: off
const elements = Dict{String,Element}([
    ["X" , "NotFound" ]              .=> Element(:X  ,"X" , "NotFound", 0,  0.00000, false, NaN);
    ["H" , "Hydrogen" ]              .=> Element(:H  ,"H" , "Hydrogen", 1,  1.00797, false, 1.1);
    ["He", "HEL", "Helium" ]         .=> Element(:He ,"He", "Helium", 2,  4.00260, false, 1.4);
    ["Li", "LI" , "Lithium" ]        .=> Element(:Li ,"Li", "Lithium", 3,  6.941, false, 1.82);
    ["Be", "BE" , "Beryllium" ]      .=> Element(:Be ,"Be", "Beryllium", 4,  9.01218, false, 1.53);
    ["B" , "Boron" ]                 .=> Element(:B  ,"B" , "Boron", 5,  10.81, false, 1.92);
    ["C" , "Carbon" ]                .=> Element(:C  ,"C" , "Carbon", 6,  12.011, false, 1.70);
    ["N" , "Nitrogen"]               .=> Element(:N  ,"N" , "Nitrogen", 7,  14.0067, false, 1.55);
    ["O" , "Oxygen" ]                .=> Element(:O  ,"O" , "Oxygen", 8,  15.9994, false, 1.52);
    ["F" , "Fluorine" ]              .=> Element(:F  ,"F" , "Fluorine", 9,  18.998403, false, 1.47);
    ["Ne", "NEO", "Neon" ]           .=> Element(:Ne ,"Ne", "Neon", 10,  20.179, false, 1.54);
    ["Na", "SOD", "NA" , "Sodium" ]  .=> Element(:Na ,"Na", "Sodium", 11,  22.98977, false, 2.27);
    ["Mg", "MG" , "Magnesium" ]      .=> Element(:Mg ,"Mg", "Magnesium", 12,  24.305, false, 1.73);
    ["Al", "AL" , "Aluminum" ]       .=> Element(:Al ,"Al", "Aluminum", 13,  26.98154, false, 1.84);
    ["Si", "SI" , "Silicon" ]        .=> Element(:Si ,"Si", "Silicon", 14,  28.0855, false, 2.10);
    ["P" , "Phosphorus" ]            .=> Element(:P  ,"P" , "Phosphorus", 15,  30.97376, false, 1.80);
    ["S" , "Sulfur" ]                .=> Element(:S  ,"S" , "Sulfur", 16,  32.06, false, 1.80);
    ["Cl", "CL" , "CLA", "Chlorine" ].=> Element(:Cl ,"Cl", "Chlorine", 17,  35.453, false, 1.75);
    ["Ar", "AR" , "Argon" ]          .=> Element(:Ar ,"Ar", "Argon", 18,  39.948, false, 1.88);
    ["K" , "POT", "Potassium" ]      .=> Element(:K  ,"K" , "Potassium", 19,  39.0983, false, 2.75);
    ["Ca", "CAL", "Calcium" ]        .=> Element(:Ca ,"Ca", "Calcium", 20,  40.08, false, 2.31);
    ["Sc", "SC" , "Scandium" ]       .=> Element(:Sc ,"Sc", "Scandium", 21,  44.9559, false, 2.11);
    ["Ti", "TI" , "Titanium" ]       .=> Element(:Ti ,"Ti", "Titanium", 22,  47.90, false, NaN);
    ["V" , "Vanadium" ]              .=> Element(:V  ,"V" , "Vanadium", 23,  50.9415, false, NaN);
    ["Cr", "CR" , "Chromium" ]       .=> Element(:Cr ,"Cr", "Chromium", 24,  51.996, false, NaN);
    ["Mn", "MN" , "Manganese" ]      .=> Element(:Mn ,"Mn", "Manganese", 25,  54.9380, false, NaN);
    ["Fe", "FE" , "Iron" ]           .=> Element(:Fe ,"Fe", "Iron", 26,  55.847, false, NaN);
    ["Co", "CO" , "Cobalt" ]         .=> Element(:Co ,"Co", "Cobalt", 27,  58.9332, false, NaN);
    ["Ni", "NI" , "Nickel" ]         .=> Element(:Ni ,"Ni", "Nickel", 28,  58.70, false, 1.63);
    ["Cu", "CU" , "Copper" ]         .=> Element(:Cu ,"Cu", "Copper", 29,  63.546, false, 1.40);
    ["Zn", "ZN" , "Zinc" ]           .=> Element(:Zn ,"Zn", "Zinc", 30,  65.38, false, 1.39);
    ["Ga", "GA" , "Gallium" ]        .=> Element(:Ga ,"Ga", "Gallium", 31,  69.72, false, 1.87);
    ["Ge", "GE" , "Germanium" ]      .=> Element(:Ge ,"Ge", "Germanium", 32,  72.59, false, 2.11);
    ["As", "AS" , "Arsenic" ]        .=> Element(:As ,"As", "Arsenic", 33,  74.9216, false, 1.85);
    ["Se", "SE" , "Selenium" ]       .=> Element(:Se ,"Se", "Selenium", 34,  78.96, false, 1.90);
    ["Br", "BR" , "Bromine" ]        .=> Element(:Br ,"Br", "Bromine", 35,  79.904, false, 1.85);
    ["Kr", "KR" , "Krypton" ]        .=> Element(:Kr ,"Kr", "Krypton", 36,  83.80, false, 2.02);
    ["Rb", "RB" , "Rubidium" ]       .=> Element(:Rb ,"Rb", "Rubidium", 37,  85.4678, false, 3.03);
    ["Sr", "SR" , "Strontium" ]      .=> Element(:Sr ,"Sr", "Strontium", 38,  87.62, false, 2.49);
    ["Y" , "Yttrium" ]               .=> Element(:Y  ,"Y" , "Yttrium", 39,  88.9059, false, NaN);
    ["Zr", "ZR" , "Zirconium" ]      .=> Element(:Zr ,"Zr", "Zirconium", 40,  91.22, false, NaN);
    ["Nb", "NB" , "Niobium" ]        .=> Element(:Nb ,"Nb", "Niobium", 41,  92.9064, false, NaN);
    ["Mo", "MO" , "Molybdenum" ]     .=> Element(:Mo ,"Mo", "Molybdenum", 42,  95.94, false, NaN);
    ["Tc", "TC" , "Technetium" ]     .=> Element(:Tc ,"Tc", "Technetium", 43,  98.0, false, NaN);
    ["Ru", "RU" , "Ruthenium" ]      .=> Element(:Ru ,"Ru", "Ruthenium", 44,  101.07, false, NaN);
    ["Rh", "RH" , "Rhodium" ]        .=> Element(:Rh ,"Rh", "Rhodium", 45,  102.9055, false, NaN);
    ["Pd", "PD" , "Palladium" ]      .=> Element(:Pd ,"Pd", "Palladium", 46,  106.4, false, 1.63);
    ["Ag", "AG" , "Silver" ]         .=> Element(:Ag ,"Ag", "Silver", 47,  107.868, false, 1.72);
    ["Cd", "CAD", "Cadmium" ]        .=> Element(:Cd ,"Cd", "Cadmium", 48,  112.41, false, 1.58);
    ["In", "IN" , "Indium" ]         .=> Element(:In ,"In", "Indium", 49,  114.82, false, 1.93);
    ["Sn", "SN" , "Tin" ]            .=> Element(:Sn ,"Sn", "Tin", 50,  118.69, false, 2.17);
    ["Sb", "SB" , "Antimony" ]       .=> Element(:Sb ,"Sb", "Antimony", 51,  121.75, false, 2.06);
    ["Te", "TE" , "Tellurium" ]      .=> Element(:Te ,"Te", "Tellurium", 52,  127.60, false, 2.06);
    ["I" , "Iodine" ]                .=> Element(:I  ,"I" , "Iodine", 53,  126.9045, false, 1.98);
    ["Xe", "XE" , "Xenon" ]          .=> Element(:Xe ,"Xe", "Xenon", 54,  131.30, false, 2.16);
    ["Cs", "CES", "Cesium" ]         .=> Element(:Cs ,"Cs", "Cesium", 55,  132.9054, false, 3.43);
    ["Ba", "BA" , "Barium" ]         .=> Element(:Ba ,"Ba", "Barium", 56,  137.33, false, 2.68);
    ["La", "LA" , "Lanthanum" ]      .=> Element(:La ,"La", "Lanthanum", 57,  138.9055, false, NaN);
    ["Ce", "CER", "Cerium" ]         .=> Element(:Ce ,"Ce", "Cerium", 58,  140.12, false, NaN);
    ["Pr", "PR" , "Praseodymium" ]   .=> Element(:Pr ,"Pr", "Praseodymium", 59,  140.9077, false, NaN);
    ["Nd", "ND" , "Neodymium" ]      .=> Element(:Nd ,"Nd", "Neodymium", 60,  144.24, false, NaN);
    ["Pm", "PM" , "Promethium" ]     .=> Element(:Pm ,"Pm", "Promethium", 61,  145, false, NaN);
    ["Sm", "SM" , "Samarium" ]       .=> Element(:Sm ,"Sm", "Samarium", 62,  150.4, false, NaN);
    ["Eu", "EU" , "Europium" ]       .=> Element(:Eu ,"Eu", "Europium", 63,  151.96, false, NaN);
    ["Gd", "GD" , "Gadolinium" ]     .=> Element(:Gd ,"Gd", "Gadolinium", 64,  157.25, false, NaN);
    ["Tb", "TB" , "Terbium" ]        .=> Element(:Tb ,"Tb", "Terbium", 65,  158.9254, false, NaN);
    ["Dy", "DY" , "Dysprosium" ]     .=> Element(:Dy ,"Dy", "Dysprosium", 66,  162.50, false, NaN);
    ["Ho", "HLM", "Holmium" ]        .=> Element(:Ho ,"Ho", "Holmium", 67,  164.9304, false, NaN);
    ["Er", "ER" , "Erbium" ]         .=> Element(:Er ,"Er", "Erbium", 68,  167.26, false, NaN);
    ["Tm", "TM" , "Thulium" ]        .=> Element(:Tm ,"Tm", "Thulium", 69,  168.9342, false, NaN);
    ["Yb", "YB" , "Ytterbium" ]      .=> Element(:Yb ,"Yb", "Ytterbium", 70,  173.04, false, NaN);
    ["Lu", "LU" , "Lutetium" ]       .=> Element(:Lu ,"Lu", "Lutetium", 71,  174.967, false, NaN);
    ["Hf", "HAF", "Hafnium" ]        .=> Element(:Hf ,"Hf", "Hafnium", 72,  178.49, false, NaN);
    ["Ta", "TA" , "Tantalum" ]       .=> Element(:Ta ,"Ta", "Tantalum", 73,  180.9479, false, NaN);
    ["W" , "W"  , "Tungsten" ]       .=> Element(:W  ,"W" , "Tungsten", 74,  183.85, false, NaN);
    ["Re", "RE" , "Rhenium" ]        .=> Element(:Re ,"Re", "Rhenium", 75,  186.207, false, NaN);
    ["Os", "OS" , "Osmium" ]         .=> Element(:Os ,"Os", "Osmium", 76,  190.2, false, NaN);
    ["Ir", "IR" , "Iridium" ]        .=> Element(:Ir ,"Ir", "Iridium", 77,  192.22, false, NaN);
    ["Pt", "PT" , "Platinum" ]       .=> Element(:Pt ,"Pt", "Platinum", 78,  195.09, false, 1.75);
    ["Au", "AU" , "Gold" ]           .=> Element(:Au ,"Au", "Gold", 79,  196.9665, false, 1.66);
    ["Hg", "MER", "Mercury" ]        .=> Element(:Hg ,"Hg", "Mercury", 80,  200.59, false, 1.55);
    ["Tl", "TL" , "Thallium" ]       .=> Element(:Tl ,"Tl", "Thallium", 81,  204.37, false, 1.96);
    ["Pb", "PB" , "Lead" ]           .=> Element(:Pb ,"Pb", "Lead", 82,  207.2, false, 2.02);
    ["Bi", "BI" , "Bismuth" ]        .=> Element(:Bi ,"Bi", "Bismuth", 83,  208.9804, false, 2.07);
    ["Po", "PO" , "Polonium" ]       .=> Element(:Po ,"Po", "Polonium", 84,  209, false, 1.97);
    ["At", "AT" , "Astatine" ]       .=> Element(:At ,"At", "Astatine", 85,  210, false, 2.02);
    ["Rn", "RN" , "Radon" ]          .=> Element(:Rn ,"Rn", "Radon", 86,  222, false, 2.20);
    ["Fr", "FR" , "Francium" ]       .=> Element(:Fr ,"Fr", "Francium", 87,  223, false, 3.48);
    ["Ra", "RA" , "Radium" ]         .=> Element(:Ra ,"Ra", "Radium", 88,  226.0254, false, 2.83);
    ["Ac", "AC" , "Actinium" ]       .=> Element(:Ac ,"Ac", "Actinium", 89,  227.0278, false, NaN);
    ["Th", "TH" , "Thorium" ]        .=> Element(:Th ,"Th", "Thorium", 90,  232.0381, false, NaN);
    ["Pa", "PA" , "Protactinium" ]   .=> Element(:Pa ,"Pa", "Protactinium", 91,  231.0359, false, NaN);
    ["U" , "U"  , "Uranium" ]        .=> Element(:U  ,"U" , "Uranium", 92,  238.029, false, 1.86);
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
PDBTools.Element(:C, InlineStrings.String3("C"), "Carbon", 6, 12.011f0, true, 1.7f0)

julia> add_element!("A2", PDBTools.elements["N"])
PDBTools.Element(:N, InlineStrings.String3("N"), "Nitrogen", 7, 14.0067f0, true, 1.55f0)

julia> element(atoms[1])
"C"

julia> element(atoms[2])
"N"

julia> mass(atoms)
26.0177f0

julia> remove_custom_elements!(); 
```

Here we repeteadly call `remove_custom_elements!()` to guarantee the proper execution of the
test codes, without any custom elements predefined.

"""
function add_element!(symbol::String, reference_element::Element; elements=elements, kargs...)
    if symbol in keys(elements)
        @warn """\n
            Element $symbol already exists. Overwriting.

        """ _file = nothing _line = nothing
    end
    properties = Dict{Symbol, Any}()
    for field in fieldnames(Element) 
        if field == :custom
            value = true
        elseif field in keys(kargs)
            value = kargs[field]
        else
            value = getfield(reference_element, field)
        end
        properties[field] = value
    end
    elements[symbol] = Element(;properties...)
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
PDBTools.Element(:N, InlineStrings.String3("N"), "Nitrogen", 7, 14.0067f0, true, 1.55f0)

julia> element(Atom(name="GN"))
"N"

julia> remove_custom_elements!();

julia> element(Atom(name="GN")) # returns `nothing`

```

Here we repeatedly call `remove_custom_elements!()` to guarantee the proper execution of the
test codes, without any custom elements predefined.

"""
function remove_custom_elements!(elements=PDBTools.elements)
    filter!(name -> !elements[name].custom, element_names)
    filter!(el -> !last(el).custom, PDBTools.elements)
    return PDBTools.elements
end

@testitem "elements" begin
    atoms = read_pdb(PDBTools.TESTPDB, "protein")
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
        @test getfield(PDBTools.elements["GN"], prop) == getfield(PDBTools.elements["N"], prop)
    end
    @test mass(Atom(name="GN")) == mass(Atom(name="N"))
    @test length(PDBTools.element_names) == 268
    remove_custom_elements!()
    @test e_ref == PDBTools.elements
    @test length(PDBTools.element_names) == 267
    @test mass(Atom(name="GN")) === nothing
    remove_custom_elements!()
end
