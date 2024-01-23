#
# List of elements with properties
#
struct Element
    symbol::Symbol
    symbol_string::String
    name::String
    atomic_number::Int
    mass::Float64
end
import Base.Broadcast.broadcastable
broadcastable(element::Element) = Ref(element)
#! format: off
const elements = Dict{String,Element}([
    ["X" , "NotFound" ]              .=> Element(:X  ,"X" , "NotFound", 0,  0.00000);
    ["H" , "Hydrogen" ]              .=> Element(:H  ,"H" , "Hydrogen", 1,  1.00797);
    ["He", "HEL", "Helium" ]         .=> Element(:He ,"He", "Helium", 2,  4.00260);
    ["Li", "LI" , "Lithium" ]        .=> Element(:Li ,"Li", "Lithium", 3,  6.941);
    ["Be", "BE" , "Beryllium" ]      .=> Element(:Be ,"Be", "Beryllium", 4,  9.01218);
    ["B" , "Boron" ]                 .=> Element(:B  ,"B" , "Boron", 5,  10.81);
    ["C" , "Carbon" ]                .=> Element(:C  ,"C" , "Carbon", 6,  12.011);
    ["N" , "Nitrogen"]               .=> Element(:N  ,"N" , "Nitrogen", 7,  14.0067);
    ["O" , "Oxygen" ]                .=> Element(:O  ,"O" , "Oxygen", 8,  15.9994);
    ["F" , "Fluorine" ]              .=> Element(:F  ,"F" , "Fluorine", 9,  18.998403);
    ["Ne", "NEO", "Neon" ]           .=> Element(:Ne ,"Ne", "Neon", 10,  20.179);
    ["Na", "SOD", "NA" , "Sodium" ]  .=> Element(:Na ,"Na", "Sodium", 11,  22.98977);
    ["Mg", "MG" , "Magnesium" ]      .=> Element(:Mg ,"Mg", "Magnesium", 12,  24.305);
    ["Al", "AL" , "Aluminum" ]       .=> Element(:Al ,"Al", "Aluminum", 13,  26.98154);
    ["Si", "SI" , "Silicon" ]        .=> Element(:Si ,"Si", "Silicon", 14,  28.0855);
    ["P" , "Phosphorus" ]            .=> Element(:P  ,"P" , "Phosphorus", 15,  30.97376);
    ["S" , "Sulfur" ]                .=> Element(:S  ,"S" , "Sulfur", 16,  32.06);
    ["Cl", "CL" , "CLA", "Chlorine" ].=> Element(:Cl ,"Cl", "Chlorine", 17,  35.453);
    ["Ar", "AR" , "Argon" ]          .=> Element(:Ar ,"Ar", "Argon", 18,  39.948);
    ["K" , "POT", "Potassium" ]      .=> Element(:K  ,"K" , "Potassium", 19,  39.0983);
    ["Ca", "CAL", "Calcium" ]        .=> Element(:Ca ,"Ca", "Calcium", 20,  40.08);
    ["Sc", "SC" , "Scandium" ]       .=> Element(:Sc ,"Sc", "Scandium", 21,  44.9559);
    ["Ti", "TI" , "Titanium" ]       .=> Element(:Ti ,"Ti", "Titanium", 22,  47.90);
    ["V" , "Vanadium" ]              .=> Element(:V  ,"V" , "Vanadium", 23,  50.9415);
    ["Cr", "CR" , "Chromium" ]       .=> Element(:Cr ,"Cr", "Chromium", 24,  51.996);
    ["Mn", "MN" , "Manganese" ]      .=> Element(:Mn ,"Mn", "Manganese", 25,  54.9380);
    ["Fe", "FE" , "Iron" ]           .=> Element(:Fe ,"Fe", "Iron", 26,  55.847);
    ["Co", "CO" , "Cobalt" ]         .=> Element(:Co ,"Co", "Cobalt", 27,  58.9332);
    ["Ni", "NI" , "Nickel" ]         .=> Element(:Ni ,"Ni", "Nickel", 28,  58.70);
    ["Cu", "CU" , "Copper" ]         .=> Element(:Cu ,"Cu", "Copper", 29,  63.546);
    ["Zn", "ZN" , "Zinc" ]           .=> Element(:Zn ,"Zn", "Zinc", 30,  65.38);
    ["Ga", "GA" , "Gallium" ]        .=> Element(:Ga ,"Ga", "Gallium", 31,  69.72);
    ["Ge", "GE" , "Germanium" ]      .=> Element(:Ge ,"Ge", "Germanium", 32,  72.59);
    ["As", "AS" , "Arsenic" ]        .=> Element(:As ,"As", "Arsenic", 33,  74.9216);
    ["Se", "SE" , "Selenium" ]       .=> Element(:Se ,"Se", "Selenium", 34,  78.96);
    ["Br", "BR" , "Bromine" ]        .=> Element(:Br ,"Br", "Bromine", 35,  79.904);
    ["Kr", "KR" , "Krypton" ]        .=> Element(:Kr ,"Kr", "Krypton", 36,  83.80);
    ["Rb", "RB" , "Rubidium" ]       .=> Element(:Rb ,"Rb", "Rubidium", 37,  85.4678);
    ["Sr", "SR" , "Strontium" ]      .=> Element(:Sr ,"Sr", "Strontium", 38,  87.62);
    ["Y" , "Yttrium" ]               .=> Element(:Y  ,"Y" , "Yttrium", 39,  88.9059);
    ["Zr", "ZR" , "Zirconium" ]      .=> Element(:Zr ,"Zr", "Zirconium", 40,  91.22);
    ["Nb", "NB" , "Niobium" ]        .=> Element(:Nb ,"Nb", "Niobium", 41,  92.9064);
    ["Mo", "MO" , "Molybdenum" ]     .=> Element(:Mo ,"Mo", "Molybdenum", 42,  95.94);
    ["Tc", "TC" , "Technetium" ]     .=> Element(:Tc ,"Tc", "Technetium", 43,  98.0);
    ["Ru", "RU" , "Ruthenium" ]      .=> Element(:Ru ,"Ru", "Ruthenium", 44,  101.07);
    ["Rh", "RH" , "Rhodium" ]        .=> Element(:Rh ,"Rh", "Rhodium", 45,  102.9055);
    ["Pd", "PD" , "Palladium" ]      .=> Element(:Pd ,"Pd", "Palladium", 46,  106.4);
    ["Ag", "AG" , "Silver" ]         .=> Element(:Ag ,"Ag", "Silver", 47,  107.868);
    ["Cd", "CAD", "Cadmium" ]        .=> Element(:Cd ,"Cd", "Cadmium", 48,  112.41);
    ["In", "IN" , "Indium" ]         .=> Element(:In ,"In", "Indium", 49,  114.82);
    ["Sn", "SN" , "Tin" ]            .=> Element(:Sn ,"Sn", "Tin", 50,  118.69);
    ["Sb", "SB" , "Antimony" ]       .=> Element(:Sb ,"Sb", "Antimony", 51,  121.75);
    ["Te", "TE" , "Tellurium" ]      .=> Element(:Te ,"Te", "Tellurium", 52,  127.60);
    ["I" , "Iodine" ]                .=> Element(:I  ,"I" , "Iodine", 53,  126.9045);
    ["Xe", "XE" , "Xenon" ]          .=> Element(:Xe ,"Xe", "Xenon", 54,  131.30);
    ["Cs", "CES", "Cesium" ]         .=> Element(:Cs ,"Cs", "Cesium", 55,  132.9054);
    ["Ba", "BA" , "Barium" ]         .=> Element(:Ba ,"Ba", "Barium", 56,  137.33);
    ["La", "LA" , "Lanthanum" ]      .=> Element(:La ,"La", "Lanthanum", 57,  138.9055);
    ["Ce", "CER", "Cerium" ]         .=> Element(:Ce ,"Ce", "Cerium", 58,  140.12);
    ["Pr", "PR" , "Praseodymium" ]   .=> Element(:Pr ,"Pr", "Praseodymium", 59,  140.9077);
    ["Nd", "ND" , "Neodymium" ]      .=> Element(:Nd ,"Nd", "Neodymium", 60,  144.24);
    ["Pm", "PM" , "Promethium" ]     .=> Element(:Pm ,"Pm", "Promethium", 61,  145);
    ["Sm", "SM" , "Samarium" ]       .=> Element(:Sm ,"Sm", "Samarium", 62,  150.4);
    ["Eu", "EU" , "Europium" ]       .=> Element(:Eu ,"Eu", "Europium", 63,  151.96);
    ["Gd", "GD" , "Gadolinium" ]     .=> Element(:Gd ,"Gd", "Gadolinium", 64,  157.25);
    ["Tb", "TB" , "Terbium" ]        .=> Element(:Tb ,"Tb", "Terbium", 65,  158.9254);
    ["Dy", "DY" , "Dysprosium" ]     .=> Element(:Dy ,"Dy", "Dysprosium", 66,  162.50);
    ["Ho", "HLM", "Holmium" ]        .=> Element(:Ho ,"Ho", "Holmium", 67,  164.9304);
    ["Er", "ER" , "Erbium" ]         .=> Element(:Er ,"Er", "Erbium", 68,  167.26);
    ["Tm", "TM" , "Thulium" ]        .=> Element(:Tm ,"Tm", "Thulium", 69,  168.9342);
    ["Yb", "YB" , "Ytterbium" ]      .=> Element(:Yb ,"Yb", "Ytterbium", 70,  173.04);
    ["Lu", "LU" , "Lutetium" ]       .=> Element(:Lu ,"Lu", "Lutetium", 71,  174.967);
    ["Hf", "HAF", "Hafnium" ]        .=> Element(:Hf ,"Hf", "Hafnium", 72,  178.49);
    ["Ta", "TA" , "Tantalum" ]       .=> Element(:Ta ,"Ta", "Tantalum", 73,  180.9479);
    ["W" , "W"  , "Tungsten" ]       .=> Element(:W  ,"W" , "Tungsten", 74,  183.85);
    ["Re", "RE" , "Rhenium" ]        .=> Element(:Re ,"Re", "Rhenium", 75,  186.207);
    ["Os", "OS" , "Osmium" ]         .=> Element(:Os ,"Os", "Osmium", 76,  190.2);
    ["Ir", "IR" , "Iridium" ]        .=> Element(:Ir ,"Ir", "Iridium", 77,  192.22);
    ["Pt", "PT" , "Platinum" ]       .=> Element(:Pt ,"Pt", "Platinum", 78,  195.09);
    ["Au", "AU" , "Gold" ]           .=> Element(:Au ,"Au", "Gold", 79,  196.9665);
    ["Hg", "MER", "Mercury" ]        .=> Element(:Hg ,"Hg", "Mercury", 80,  200.59);
    ["Tl", "TL" , "Thallium" ]       .=> Element(:Tl ,"Tl", "Thallium", 81,  204.37);
    ["Pb", "PB" , "Lead" ]           .=> Element(:Pb ,"Pb", "Lead", 82,  207.2);
    ["Bi", "BI" , "Bismuth" ]        .=> Element(:Bi ,"Bi", "Bismuth", 83,  208.9804);
    ["Po", "PO" , "Polonium" ]       .=> Element(:Po ,"Po", "Polonium", 84,  209);
    ["At", "AT" , "Astatine" ]       .=> Element(:At ,"At", "Astatine", 85,  210);
    ["Rn", "RN" , "Radon" ]          .=> Element(:Rn ,"Rn", "Radon", 86,  222);
    ["Fr", "FR" , "Francium" ]       .=> Element(:Fr ,"Fr", "Francium", 87,  223);
    ["Ra", "RA" , "Radium" ]         .=> Element(:Ra ,"Ra", "Radium", 88,  226.0254);
    ["Ac", "AC" , "Actinium" ]       .=> Element(:Ac ,"Ac", "Actinium", 89,  227.0278);
    ["Th", "TH" , "Thorium" ]        .=> Element(:Th ,"Th", "Thorium", 90,  232.0381);
    ["Pa", "PA" , "Protactinium" ]   .=> Element(:Pa ,"Pa", "Protactinium", 91,  231.0359);
    ["U" , "U"  , "Uranium" ]        .=> Element(:U  ,"U" , "Uranium", 92,  238.029);
])
#! format: on

# Sort element_names by name for faster searching
const element_names = sort(collect(keys(elements)))

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
end
