#
# Ramachandran plots: scatter(::Ramachandran)
#

"""
    scatter(ram::Ramachandran; kwargs...)

Ramachandran plot of the ϕ and ψ angles of a [`Ramachandran`](@ref) object. Returns a
`Makie.Figure`.

!!! note
    This function requires a Makie backend, e.g. `using CairoMakie`.

# Keyword arguments
- `xlabel = "ϕ (degrees)"`, `ylabel = "ψ (degrees)"`, `title = ""`: labels of the axes, and title.
- `size = (450, 450)`, `font = nothing`, `fontsize = 14`: size of the figure, font
  (default: New Computer Modern, distributed with Makie) and font size.
- Other keywords are passed to `Makie.scatter!` (e.g. `color`, `markersize`).

# Example

```julia-repl
julia> using PDBTools, LazyArtifacts, CairoMakie

julia> prot = read_pdb(PDBTools.test_file("structure.pdb"), "protein");

julia> scatter(Ramachandran(prot))
```
"""
function Makie.scatter(ram::Ramachandran; size = (450, 450), font = nothing, fontsize = 14, kwargs...)
    fig = _figure(; size, font, fontsize)
    Makie.scatter(fig[1, 1], ram; kwargs...)
    return fig
end

function Makie.scatter(gp::GridPos, ram::Ramachandran;
    xlabel::AbstractString = "ϕ (degrees)",
    ylabel::AbstractString = "ψ (degrees)",
    title::AbstractString = "",
    kwargs...,
)
    gl = GridLayout(gp)
    ax = Axis(gl[1, 1]; xlabel, ylabel, title, aspect = DataAspect(), xticks = -180:45:180, yticks = -180:45:180)
    Makie.scatter!(ax, ram; kwargs...)
    limits!(ax, -180, 180, -180, 180)
    return gl
end

Makie.scatter!(ax::Axis, ram::Ramachandran; kwargs...) = Makie.scatter!(ax, ram.phi, ram.psi; kwargs...)

@testitem "scatter(::Ramachandran)" begin
    using PDBTools, LazyArtifacts
    import CairoMakie
    prot = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    ram = Ramachandran(prot)
    fig = CairoMakie.scatter(ram; color = :red, markersize = 6)
    @test fig isa CairoMakie.Figure
    tmpplot = tempname() * ".png"
    CairoMakie.save(tmpplot, fig)
    @test isfile(tmpplot)
    fig = CairoMakie.Figure()
    @test CairoMakie.scatter(fig[1, 1], ram) isa CairoMakie.GridLayout
    @test CairoMakie.scatter!(CairoMakie.Axis(fig[1, 2]), ram) isa CairoMakie.AbstractPlot
end
