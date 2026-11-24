#
# Deprecated: plots with Plots.jl. The plots of PDBTools are made with Makie (see the
# MakiePlotting extension). These methods are kept, without new features, for
# compatibility, and will be removed in a future breaking release.
#
module PlotsPlotting
    using TestItems: @testitem
    using PDBTools
    import Plots
    include("./contacts.jl")
    include("./ramachandran.jl")
end
