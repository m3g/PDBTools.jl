function PDBTools.plot_hammeraitoff(
    atoms::AbstractVector{<:PDBTools.Atom},
    values;
    nbins::Integer=100,
    statistic::Symbol=:mean,
    colormap=:tempo,
    clims=nothing,
    contours::Symbol=:none,
    groupfun=identity,
    contour_color=:grey30,
    colorbar_title::AbstractString="",
    title::AbstractString="",
    outfile=nothing,
)
    if !(contours in (:none, :chain, :chains, :group))
        throw(ArgumentError("contours must be :none, :chain (or :chains) or :group, got $contours"))
    end
    length(values) == length(atoms) ||
        throw(ArgumentError("values has $(length(values)) elements, but there are $(length(atoms)) atoms"))

    # Project once and reuse the result for the heat map and the contours
    x_hammer, y_hammer = PDBTools.hammer(atoms)
    xbin, ybin, cbins = PDBTools.hand_made_histogram(x_hammer, y_hammer, values; nbins, statistic)

    p = Plots.contourf(
        xbin, ybin, cbins;
        color=colormap,
        clims=clims,
        xticks=[],
        yticks=[],
        lw=0.0,
        label="",
        colorbar_title=colorbar_title,
        title=title,
        framestyle=:box,
        fontfamily="Computer Modern",
    )

    if contours != :none
        # Label of each atom: its chain, or the group of its chain
        labels = if contours == :group
            [string(groupfun(PDBTools.chain(at))) for at in atoms]
        else
            [string(PDBTools.chain(at)) for at in atoms]
        end

        # The color scale is shared by all series of the plot. To keep the
        # contour masks from changing it, the masks use only the lowest and
        # highest values of the color scale, and `clims` is passed to every
        # `contour!` call.
        lo, hi = isnothing(clims) ? extrema(cbins) : clims

        # One contour per label: the mask is `hi` for atoms with that label
        # and `lo` for all the others, and the contour is drawn at the border.
        for label in unique(labels)
            mask_values = [l == label ? hi : lo for l in labels]
            xb, yb, mask = PDBTools.hand_made_histogram(x_hammer, y_hammer, mask_values; nbins, statistic=:mean)
            Plots.contour!(
                p, xb, yb, mask;
                levels=1,
                clims=(lo, hi),
                color=contour_color,
                lw=0.5,
                label="",
            )
        end
    end

    isnothing(outfile) || Plots.savefig(p, outfile)
    return p
end

@testitem "plot_hammeraitoff" begin
    using PDBTools, LazyArtifacts
    import Plots
    ENV["GKSwstype"] = "nul"
    atoms = read_pdb(PDBTools.test_file("6co8.pdb"), "protein")
    n = length(atoms)
    for contours in (:none, :chain)
        tmpfile = tempname() * ".png"
        p = plot_hammeraitoff(atoms, beta.(atoms); nbins = 30, contours, outfile = tmpfile)
        @test p isa Plots.Plot
        @test isfile(tmpfile)
    end
    p = plot_hammeraitoff(atoms, ones(n); nbins = 30, contours = :group, groupfun = c -> c == "A" ? "A" : "other")
    @test p isa Plots.Plot
    @test_throws ArgumentError plot_hammeraitoff(atoms, ones(n); contours = :residue)
    @test_throws ArgumentError plot_hammeraitoff(atoms, ones(n - 1))
end
