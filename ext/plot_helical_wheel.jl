function PDBTools.helical_wheel(
    sequence::AbstractString,
    values::AbstractVector{<:Real};
    residue_numbers::AbstractVector{<:Integer} = 1:length(sequence),
    normalize::Bool = true,
    colormap = :tempo,
    colorbar_title::AbstractString = "",
    clims = nothing,
    hydrophobicity::AbstractDict = PDBTools.FAUCHERE_PLISKA,
    show_moment::Bool = true,
)
    seq = collect(uppercase(sequence))
    n = length(seq)
    length(values) == n ||
        throw(ArgumentError("values has $(length(values)) elements, but the sequence has $n residues."))
    length(residue_numbers) == n ||
        throw(ArgumentError("residue_numbers has $(length(residue_numbers)) elements, but the sequence has $n residues."))

    # Mean hydrophobicity and hydrophobic moment (also validates the sequence)
    moment = PDBTools.hydrophobic_moment(sequence; hydrophobicity)
    μx, μy = moment.μx, moment.μy

    # Optional min-max normalization
    weights = float.(collect(values))
    if normalize
        min_w, max_w = extrema(weights)
        weights = (weights .- min_w) ./ (max_w - min_w + 1e-9)
    end

    # Geometry: 100° per residue; new ring every 18 residues (5 turns)
    turn_residues = 18
    radius_step = 0.25
    angles = PDBTools._helical_wheel_angles(n)
    radii = [1 + div(i, turn_residues) * radius_step for i in 0:n-1]

    # Rotate the diagram so that the moment points down
    rotation_angle = atan(μy, μx) + π / 2
    rotated_angles = angles .- rotation_angle
    μx_rot = μx * cos(rotation_angle) + μy * sin(rotation_angle)
    μy_rot = -μx * sin(rotation_angle) + μy * cos(rotation_angle)

    x = cos.(rotated_angles) .* radii
    y = sin.(rotated_angles) .* radii

    p = Plots.plot(x, y; lw = 1, lc = :gray10, alpha = 0.6, label = "", fontfamily = "Computer Modern")
    Plots.scatter!(p, x, y;
        marker_z = weights,
        markersize = 14,
        color = Plots.cgrad(colormap),
        colorbar = :left,
        colorbar_title = colorbar_title,
        label = "",
    )
    isnothing(clims) || Plots.plot!(p; clims = clims)

    for i in 1:n
        Plots.annotate!(p, x[i], y[i], Plots.text(seq[i], :black, :center))
    end
    Plots.annotate!(p, x[1], y[1] - 0.16, Plots.text("$(first(residue_numbers))", :gray10, :center, 4))
    Plots.annotate!(p, x[end], y[end] - 0.16, Plots.text("$(last(residue_numbers))", :gray10, :center, 4))

    if show_moment
        Plots.quiver!(p, [0], [0]; quiver = ([μx_rot], [μy_rot]), color = :gray30, linewidth = 1, label = "")
    end

    lim = maximum(radii) + 0.5
    Plots.plot!(p;
        xlim = (-lim, lim), ylim = (-lim, lim),
        aspect_ratio = 1, grid = false, legend = false, frame = :box, ticks = [],
        background_color = :transparent, dpi = 300,
    )

    return p, (mean_hydro = moment.mean_hydro, μ = moment.μ)
end

function PDBTools.helical_wheel(
    atoms::AbstractVector{<:PDBTools.Atom};
    aggregate = PDBTools._mean,
    property = PDBTools.beta,
    kwargs...,
)
    values, sequence, numbers = PDBTools.residue_values(atoms; aggregate, property)
    return PDBTools.helical_wheel(sequence, values; residue_numbers = numbers, kwargs...)
end

@testitem "helical_wheel" begin
    using PDBTools, LazyArtifacts
    import Plots
    ENV["GKSwstype"] = "nul"

    seq = "LKKLLKLLKKLLKL"
    n = length(seq)
    p, info = helical_wheel(seq, collect(1.0:n); colorbar_title = "CN")
    @test p isa Plots.Plot
    @test info.μ ≈ hydrophobic_moment(seq).μ
    @test info.mean_hydro ≈ hydrophobic_moment(seq).mean_hydro
    tmpfile = tempname() * ".png"
    Plots.savefig(p, tmpfile)
    @test isfile(tmpfile)
    # more than 18 residues (second ring), fixed color limits, no arrow
    p, _ = helical_wheel(repeat(seq, 2), rand(2n); normalize = false, clims = (0, 1), show_moment = false)
    @test p isa Plots.Plot
    @test_throws ArgumentError helical_wheel(seq, ones(n - 1))
    @test_throws ArgumentError helical_wheel(seq, ones(n); residue_numbers = 1:n-1)
    @test_throws ArgumentError helical_wheel("LKXL", ones(4))

    # from a selection of atoms
    pdb = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    ch = chain(first(pdb))
    atoms = select(pdb, at -> chain(at) == ch && resnum(at) <= resnum(first(pdb)) + 17)
    p, info = helical_wheel(atoms; aggregate = sum)
    @test p isa Plots.Plot
end
