#
# Plots of PDBTools, with Makie (any backend: CairoMakie, GLMakie, WGLMakie).
#
# Each plotted type has three methods (for `f` = `plot`, `heatmap` or `scatter`):
#
# - `f(x; kwargs...)`: the figure, with its color bar, returned as a `Makie.Figure`. This
#   is the method most users need.
# - `f(fig[i, j], x; kwargs...)`: draws the plot (and its color bar) into a position of an
#   existing figure, for composite figures. Returns the `GridLayout` with the plot.
# - `f!(ax, x; kwargs...)`: draws the plot into an existing `Axis` (without color bar).
#   Returns the plot object. Not available for `SecondaryStructureDiagram`, which has one axis per row.
#
module MakiePlotting

using TestItems: @testitem
using PDBTools
using PDBTools: Atom, SSData, SecondaryStructureDiagram, HelicalWheel, HammerAitoffMap, ContactMap, Ramachandran
import Makie
using SparseArrays: findnz
using Makie: Figure, Axis, Colorbar, GridLayout, GridPosition, GridSubposition, Point2f, Vec2f,
    poly!, lines!, scatter!, text!, contour!, contourf!, arrows2d!,
    limits!, hidespines!, hidedecorations!, rowgap!, resize_to_layout!, cgrad, DataAspect

const GridPos = Union{GridPosition,GridSubposition}

# Fonts of the figures: New Computer Modern, distributed with Makie, unless a font is given
function _fonts(font)
    isnothing(font) || return (; regular = font, bold = font)
    return (;
        regular = Makie.texfont(:regular), bold = Makie.texfont(:bold),
        italic = Makie.texfont(:italic), bold_italic = Makie.texfont(:bolditalic),
    )
end

function _check_size(size)
    length(size) == 2 && all(>(0), size) ||
        throw(ArgumentError("size must be a tuple (width, height) of positive numbers, got $size"))
    return Float64.(Tuple(size))
end

# Figure of the `f(x; kwargs...)` methods
_figure(; size, font = nothing, fontsize = 14, kwargs...) =
    Figure(; size = _check_size(size), fonts = _fonts(font), fontsize, kwargs...)

include("./secondary_structure_diagram.jl")
include("./hammer_aitoff.jl")
include("./helical_wheel.jl")
include("./contacts.jl")
include("./ramachandran.jl")

end # module MakiePlotting
