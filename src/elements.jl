#
# List of elements with properties
#
struct Element
    symbol::String
    name::String
    atomic_number::Int
    mass::Float64
end
import Base.Broadcast.broadcastable
broadcastable(element::Element) = Ref(element)
const elements = Dict{String,Element}([
    ["X" , "NotFound" ]              .=> Element("X" , "NotFound", 0,  0.00000);
    ["H" , "Hydrogen" ]              .=> Element("H" , "Hydrogen", 1,  1.00797);
    ["He", "HEL", "Helium" ]         .=> Element("He", "Helium", 2,  4.00260);
    ["Li", "LI" , "Lithium" ]        .=> Element("Li", "Lithium", 3,  6.941);
    ["Be", "BE" , "Beryllium" ]      .=> Element("Be", "Beryllium", 4,  9.01218);
    ["B" , "Boron" ]                 .=> Element("B" , "Boron", 5,  10.81);
    ["C" , "Carbon" ]                .=> Element("C" , "Carbon", 6,  12.011);
    ["N" , "Nitrogen"]               .=> Element("N" , "Nitrogen", 7,  14.0067);
    ["O" , "Oxygen" ]                .=> Element("O" , "Oxygen", 8,  15.9994);
    ["F" , "Fluorine" ]              .=> Element("F" , "Fluorine", 9,  18.998403);
    ["Ne", "NEO", "Neon" ]           .=> Element("Ne", "Neon", 10,  20.179);
    ["Na", "SOD", "NA" , "Sodium" ]  .=> Element("Na", "Sodium", 11,  22.98977);
    ["Mg", "MG" , "Magnesium" ]      .=> Element("Mg", "Magnesium", 12,  24.305);
    ["Al", "AL" , "Aluminum" ]       .=> Element("Al", "Aluminum", 13,  26.98154);
    ["Si", "SI" , "Silicon" ]        .=> Element("Si", "Silicon", 14,  28.0855);
    ["P" , "Phosphorus" ]            .=> Element("P" , "Phosphorus", 15,  30.97376);
    ["S" , "Sulfur" ]                .=> Element("S" , "Sulfur", 16,  32.06);
    ["Cl", "CL" , "CLA", "Chlorine" ].=> Element("Cl", "Chlorine", 17,  35.453);
    ["Ar", "AR" , "Argon" ]          .=> Element("Ar", "Argon", 18,  39.948);
    ["K" , "POT", "Potassium" ]      .=> Element("K" , "Potassium", 19,  39.0983);
    ["Ca", "CAL", "Calcium" ]        .=> Element("Ca", "Calcium", 20,  40.08);
    ["Sc", "SC" , "Scandium" ]       .=> Element("Sc", "Scandium", 21,  44.9559);
    ["Ti", "TI" , "Titanium" ]       .=> Element("Ti", "Titanium", 22,  47.90);
    ["V" , "Vanadium" ]              .=> Element("V" , "Vanadium", 23,  50.9415);
    ["Cr", "CR" , "Chromium" ]       .=> Element("Cr", "Chromium", 24,  51.996);
    ["Mn", "MN" , "Manganese" ]      .=> Element("Mn", "Manganese", 25,  54.9380);
    ["Fe", "FE" , "Iron" ]           .=> Element("Fe", "Iron", 26,  55.847);
    ["Co", "CO" , "Cobalt" ]         .=> Element("Co", "Cobalt", 27,  58.9332);
    ["Ni", "NI" , "Nickel" ]         .=> Element("Ni", "Nickel", 28,  58.70);
    ["Cu", "CU" , "Copper" ]         .=> Element("Cu", "Copper", 29,  63.546);
    ["Zn", "ZN" , "Zinc" ]           .=> Element("Zn", "Zinc", 30,  65.38);
    ["Ga", "GA" , "Gallium" ]        .=> Element("Ga", "Gallium", 31,  69.72);
    ["Ge", "GE" , "Germanium" ]      .=> Element("Ge", "Germanium", 32,  72.59);
    ["As", "AS" , "Arsenic" ]        .=> Element("As", "Arsenic", 33,  74.9216);
    ["Se", "SE" , "Selenium" ]       .=> Element("Se", "Selenium", 34,  78.96);
    ["Br", "BR" , "Bromine" ]        .=> Element("Br", "Bromine", 35,  79.904);
    ["Kr", "KR" , "Krypton" ]        .=> Element("Kr", "Krypton", 36,  83.80);
    ["Rb", "RB" , "Rubidium" ]       .=> Element("Rb", "Rubidium", 37,  85.4678);
    ["Sr", "SR" , "Strontium" ]      .=> Element("Sr", "Strontium", 38,  87.62);
    ["Y" , "Yttrium" ]               .=> Element("Y" , "Yttrium", 39,  88.9059);
    ["Zr", "ZR" , "Zirconium" ]      .=> Element("Zr", "Zirconium", 40,  91.22);
    ["Nb", "NB" , "Niobium" ]        .=> Element("Nb", "Niobium", 41,  92.9064);
    ["Mo", "MO" , "Molybdenum" ]     .=> Element("Mo", "Molybdenum", 42,  95.94);
    ["Tc", "TC" , "Technetium" ]     .=> Element("Tc", "Technetium", 43,  98.0);
    ["Ru", "RU" , "Ruthenium" ]      .=> Element("Ru", "Ruthenium", 44,  101.07);
    ["Rh", "RH" , "Rhodium" ]        .=> Element("Rh", "Rhodium", 45,  102.9055);
    ["Pd", "PD" , "Palladium" ]      .=> Element("Pd", "Palladium", 46,  106.4);
    ["Ag", "AG" , "Silver" ]         .=> Element("Ag", "Silver", 47,  107.868);
    ["Cd", "CAD", "Cadmium" ]        .=> Element("Cd", "Cadmium", 48,  112.41);
    ["In", "IN" , "Indium" ]         .=> Element("In", "Indium", 49,  114.82);
    ["Sn", "SN" , "Tin" ]            .=> Element("Sn", "Tin", 50,  118.69);
    ["Sb", "SB" , "Antimony" ]       .=> Element("Sb", "Antimony", 51,  121.75);
    ["Te", "TE" , "Tellurium" ]      .=> Element("Te", "Tellurium", 52,  127.60);
    ["I" , "Iodine" ]                .=> Element("I" , "Iodine", 53,  126.9045);
    ["Xe", "XE" , "Xenon" ]          .=> Element("Xe", "Xenon", 54,  131.30);
    ["Cs", "CES", "Cesium" ]         .=> Element("Cs", "Cesium", 55,  132.9054);
    ["Ba", "BA" , "Barium" ]         .=> Element("Ba", "Barium", 56,  137.33);
    ["La", "LA" , "Lanthanum" ]      .=> Element("La", "Lanthanum", 57,  138.9055);
    ["Ce", "CER", "Cerium" ]         .=> Element("Ce", "Cerium", 58,  140.12);
    ["Pr", "PR" , "Praseodymium" ]   .=> Element("Pr", "Praseodymium", 59,  140.9077);
    ["Nd", "ND" , "Neodymium" ]      .=> Element("Nd", "Neodymium", 60,  144.24);
    ["Pm", "PM" , "Promethium" ]     .=> Element("Pm", "Promethium", 61,  145);
    ["Sm", "SM" , "Samarium" ]       .=> Element("Sm", "Samarium", 62,  150.4);
    ["Eu", "EU" , "Europium" ]       .=> Element("Eu", "Europium", 63,  151.96);
    ["Gd", "GD" , "Gadolinium" ]     .=> Element("Gd", "Gadolinium", 64,  157.25);
    ["Tb", "TB" , "Terbium" ]        .=> Element("Tb", "Terbium", 65,  158.9254);
    ["Dy", "DY" , "Dysprosium" ]     .=> Element("Dy", "Dysprosium", 66,  162.50);
    ["Ho", "HLM", "Holmium" ]        .=> Element("Ho", "Holmium", 67,  164.9304);
    ["Er", "ER" , "Erbium" ]         .=> Element("Er", "Erbium", 68,  167.26);
    ["Tm", "TM" , "Thulium" ]        .=> Element("Tm", "Thulium", 69,  168.9342);
    ["Yb", "YB" , "Ytterbium" ]      .=> Element("Yb", "Ytterbium", 70,  173.04);
    ["Lu", "LU" , "Lutetium" ]       .=> Element("Lu", "Lutetium", 71,  174.967);
    ["Hf", "HAF", "Hafnium" ]        .=> Element("Hf", "Hafnium", 72,  178.49);
    ["Ta", "TA" , "Tantalum" ]       .=> Element("Ta", "Tantalum", 73,  180.9479);
    ["W" , "W"  , "Tungsten" ]       .=> Element("W" , "Tungsten", 74,  183.85);
    ["Re", "RE" , "Rhenium" ]        .=> Element("Re", "Rhenium", 75,  186.207);
    ["Os", "OS" , "Osmium" ]         .=> Element("Os", "Osmium", 76,  190.2);
    ["Ir", "IR" , "Iridium" ]        .=> Element("Ir", "Iridium", 77,  192.22);
    ["Pt", "PT" , "Platinum" ]       .=> Element("Pt", "Platinum", 78,  195.09);
    ["Au", "AU" , "Gold" ]           .=> Element("Au", "Gold", 79,  196.9665);
    ["Hg", "MER", "Mercury" ]        .=> Element("Hg", "Mercury", 80,  200.59);
    ["Tl", "TL" , "Thallium" ]       .=> Element("Tl", "Thallium", 81,  204.37);
    ["Pb", "PB" , "Lead" ]           .=> Element("Pb", "Lead", 82,  207.2);
    ["Bi", "BI" , "Bismuth" ]        .=> Element("Bi", "Bismuth", 83,  208.9804);
    ["Po", "PO" , "Polonium" ]       .=> Element("Po", "Polonium", 84,  209);
    ["At", "AT" , "Astatine" ]       .=> Element("At", "Astatine", 85,  210);
    ["Rn", "RN" , "Radon" ]          .=> Element("Rn", "Radon", 86,  222);
    ["Fr", "FR" , "Francium" ]       .=> Element("Fr", "Francium", 87,  223);
    ["Ra", "RA" , "Radium" ]         .=> Element("Ra", "Radium", 88,  226.0254);
    ["Ac", "AC" , "Actinium" ]       .=> Element("Ac", "Actinium", 89,  227.0278);
    ["Th", "TH" , "Thorium" ]        .=> Element("Th", "Thorium", 90,  232.0381);
    ["Pa", "PA" , "Protactinium" ]   .=> Element("Pa", "Protactinium", 91,  231.0359);
    ["U" , "U"  , "Uranium" ]        .=> Element("U" , "Uranium", 92,  238.029);
])


# Function that tries to match an atom name, PDB style, with one of the
# element names
const element_names = sort(collect(keys(elements)))
function match_element(name::String)
    # if there is match, just return the name
    iel = searchsortedfirst(element_names, name)
    if iel <= length(element_names) && name == element_names[iel]
        return name
    end
    # Check if the first character is number
    i0 = 1 + isdigit(first(name))    
    imatch = searchsortedfirst(element_names, name[i0:i0]; by=x -> x[1])
    lmatch = searchsortedlast(element_names, name[i0:i0]; by=x -> x[1])
    for iel in imatch:lmatch
        el = element_names[iel]
        if el == name[i0:i0+length(el)-1]
            return el
        end
    end
    return nothing
end

#
# Retrive index of element in elements list from name. 
#
atomic_number(name::String) = isnothing(match_element(name)) ? nothing : elements[match_element(name)].atomic_number
mass(name::String) = isnothing(match_element(name)) ? nothing : elements[match_element(name)].mass
element(name::String) = isnothing(match_element(name)) ? nothing : elements[match_element(name)].symbol
element_name(name::String) = isnothing(match_element(name)) ? nothing : elements[match_element(name)].name

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
#    @test mass.(select(atoms, "residue = 1")) == 
#voltar
end
