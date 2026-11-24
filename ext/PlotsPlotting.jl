module PlotsPlotting
    using TestItems: @testitem
    using PDBTools
    import Plots
    include("./plot_contacts.jl")
    include("./plot_ramachandran.jl")
    include("./plot_hammer_aitoff.jl")
    include("./plot_helical_wheel.jl")
end
