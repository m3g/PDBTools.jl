import Plots: scatter

"""
    scatter(ram::Ramachandran; kargs...)

Creates a Ramachandran plot given a `Ramachandran` object.

# Arguments

- `ram::Ramachandran`: the Ramachandran object, containing ϕ and ψ angles, 
  resulting from the the `Ramachandan` function.

All other arguments are default keywords of `Plots.scatter` and can be adjusted to customize the plot.

# Example 

```jldoctest
julia> using PDBTools, Plots

julia> prot = read_pdb(PDBTools.TESTPDB, "protein");

julia> ram = Ramachandran(prot)
Ramachandran data: phi, psi vectors with 102 angles.

julia> # plt = scatter(map) # uncomment to plot
```

"""
function scatter(ram::Ramachandran;
    xlabel="Φ (degrees)",
    ylabel="Ψ (degrees)",
    ratio=1,
    xlims=(-180, 180),
    ylims=(-180, 180),
    framestyle=:box,
    label=:none,
    size=(400,400),
    xticks=(-180:45:180),
    yticks=(-180:45:180),
    fontfamily="Computer Modern",
    kargs...
)
    Plots.scatter(ram.phi, ram.psi;
        xlabel,
        ylabel,
        ratio,
        xlims,
        ylims,
        framestyle,
        label,
        size,
        xticks,
        yticks,
        fontfamily, 
        kargs...
    )
end

@testitem "ramachandran plots" begin
    using PDBTools
    using Plots
    tmpplot = tempname()*".png"
    prot = read_pdb(PDBTools.TESTPDB, "protein")
    ram = Ramachandran(prot)
    plt = scatter(ram)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
end
