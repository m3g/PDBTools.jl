include("./aa_properties.jl")
include("./which_natural_aminoacid.jl")
include("./selmacros.jl")
include("./select.jl")

export isacidic, isaliphatic, isaromatic, isbackbone, isbasic,
       ischarged, ishydrophobic, isneutral, isnonpolar, ispolar,
       isprotein, issidechain, iswater

export select, selindex
