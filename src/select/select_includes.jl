include("./aa_properties.jl")
include("./which_natural_aminoacid.jl")
include("./isacidic.jl")
include("./isaliphatic.jl")
include("./isaromatic.jl")
include("./isbackbone.jl")
include("./isbasic.jl")
include("./ischarged.jl")
include("./ishydrophobic.jl")
include("./isneutral.jl")
include("./isnonpolar.jl")
include("./ispolar.jl")
include("./isprotein.jl")
include("./issidechain.jl")
include("./iswater.jl")
include("./select.jl")

export isacidic, isaliphatic, isaromatic, isbackbone, isbasic,
       ischarged, ishydrophobic, isneutral, isnonpolar, ispolar,
       isprotein, issidechain, iswater

export select, selindex
