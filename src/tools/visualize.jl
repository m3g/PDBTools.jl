#
# Interactive visualization of structures using 3Dmol.js (https://3dmol.csb.pitt.edu)
#
const _3DMOL_URL = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.5/3Dmol-min.js"

"""
    StructureView

Object returned by [`visualize`](@ref). It is displayed as an interactive 3D view in
environments that render HTML (VS Code plot pane, Pluto, Jupyter, Documenter pages),
and can be saved as a standalone HTML file with `save(filename, view)`.

"""
struct StructureView
    html::String
    natoms::Int
    width::String
    height::String
end

"""
    visualize(atoms::AbstractVector{<:Atom}, [selection]; kargs...)

Returns an interactive 3D view of the structure, rendered with [3Dmol.js](https://3dmol.csb.pitt.edu).
The view is displayed in environments that render HTML (VS Code plot pane, Pluto, Jupyter notebooks,
Documenter pages). It can be saved as a standalone HTML file with `save("view.html", view)`.

Secondary structures for the cartoon representation are computed by 3Dmol.js.
The optional `selection` (a string or function) restricts the atoms shown.

The atoms are classified as polymer (protein or nucleic acid), water, or other (ligands, ions, etc.).
The main `style` applies to polymer atoms, and to ligands if `style` is not `:cartoon`.

# Keyword arguments

- `style::Symbol=:cartoon`: `:cartoon`, `:sticks`, `:ballandstick`, `:spheres`, or `:lines`.
- `color=:auto`: one of `:chain`, `:ss` (secondary structure), `:spectrum` (cartoon only), `:element`,
  `:residue`, or a color string (e.g. `"red"` or `"#ff0000"`). With `:auto`, cartoons are colored by
  chain (or by spectrum if there is a single chain), and other styles by element.
- `color_by::AbstractVector{<:Real}=nothing`: a value per atom of `atoms` (before selection),
  used to color the atoms with a gradient. Overrides `color`.
- `colormap::Symbol=:rwb`: gradient used with `color_by`: `:rwb` (red-white-blue), `:roygb`, or `:sinebow`.
- `color_range=nothing`: `(min, max)` limits of the gradient. Defaults to the range of `color_by`.
- `ligands::Bool=true`: show non-polymer, non-water atoms (as sticks if `style == :cartoon`).
- `water=:auto`: show water molecules (`true` or `false`). With `:auto`, water is shown only if there
  is nothing else to show (no polymer or other molecules).
- `surface::Bool=false`: add a transparent molecular surface for the polymer atoms.
- `hover::Bool=true`: show atom labels when the mouse is over the atoms.
- `unitcell=nothing`: draw the edges of the periodic box. A 3x3 matrix with the lattice vectors 
  as columns (as returned by [`read_unitcell`](@ref)), or a vector with the box sides, for orthorhombic boxes.
- `unitcell_origin=(0, 0, 0)`: position of the origin (corner) of the box, or `:center` to center the box 
  on the geometric center of the atoms.
- `unitcell_color="gray"`: color of the box edges.
- `background="white"`: background color.
- `width="100%"`, `height=400`: size of the view, in pixels if given as numbers.

# Example

```julia
using PDBTools
atoms = read_pdb(PDBTools.TESTPDB, "protein")
visualize(atoms)
visualize(atoms; style=:sticks, color=:residue)
visualize(atoms; color_by=beta.(atoms), colormap=:roygb)
v = visualize(atoms, "chain A"; surface=true)
visualize(read_pdb(PDBTools.TESTPBC); unitcell=read_unitcell(PDBTools.TESTPBC))
save("view.html", v)
```

"""
function visualize(
    atoms::AbstractVector{<:Atom},
    selection::AbstractString;
    kargs...
)
    visualize(atoms, parse_query(selection); kargs...)
end

function visualize(
    atoms::AbstractVector{<:Atom},
    selection_function::Function=all;
    style::Symbol=:cartoon,
    color::Union{Symbol,AbstractString}=:auto,
    color_by::Union{Nothing,AbstractVector{<:Real}}=nothing,
    colormap::Symbol=:rwb,
    color_range::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    ligands::Bool=true,
    water::Union{Bool,Symbol}=:auto,
    surface::Bool=false,
    hover::Bool=true,
    unitcell::Union{Nothing,AbstractVector{<:Real},AbstractMatrix{<:Real}}=nothing,
    unitcell_origin::Union{Symbol,AbstractVector{<:Real},Tuple{Vararg{Real,3}}}=(0, 0, 0),
    unitcell_color::AbstractString="gray",
    background::AbstractString="white",
    width::Union{Integer,AbstractString}="100%",
    height::Union{Integer,AbstractString}=400,
)
    style in (:cartoon, :sticks, :ballandstick, :spheres, :lines) ||
        throw(ArgumentError("style must be one of :cartoon, :sticks, :ballandstick, :spheres, or :lines. Got: :$style"))
    water isa Bool || water == :auto ||
        throw(ArgumentError("water must be true, false, or :auto. Got: :$water"))
    if !isnothing(unitcell) && !(size(unitcell) in ((3,), (3, 3)))
        throw(ArgumentError("unitcell must be a 3x3 matrix or a vector of length 3. Got size: $(size(unitcell))"))
    end
    unitcell_origin isa Symbol && unitcell_origin != :center &&
        throw(ArgumentError("unitcell_origin must be a vector of length 3 or :center. Got: :$unitcell_origin"))
    unitcell_origin isa AbstractVector && length(unitcell_origin) != 3 &&
        throw(ArgumentError("unitcell_origin must be a vector of length 3 or :center. Got length: $(length(unitcell_origin))"))
    colormap in (:rwb, :roygb, :sinebow) ||
        throw(ArgumentError("colormap must be one of :rwb, :roygb, or :sinebow. Got: :$colormap"))
    if !isnothing(color_by) && length(color_by) != length(atoms)
        throw(ArgumentError("color_by must have the same length as atoms ($(length(atoms))). Got: $(length(color_by))"))
    end
    # mmCIF text of the selected atoms, and group of each atom (0: polymer, 1: other, 2: water), which is
    # used in the 3Dmol.js selections through the atom index (atoms are kept in the same order)
    cif = IOBuffer()
    println(cif, "data_PDBTools")
    println(cif, "loop_")
    for field in ("group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id", "auth_asym_id",
        "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv")
        println(cif, "_atom_site.", field)
    end
    groups = IOBuffer()
    chains = Set{String}()
    natoms = 0
    npolymer = 0
    nbackbone = 0 # CA and P atoms, required for the cartoon representation
    nother = 0
    vmin, vmax = Inf, -Inf
    geometric_center = zeros(3)
    for (i, atom) in enumerate(atoms)
        selection_function(atom) || continue
        natoms += 1
        geometric_center .+= (atom.x, atom.y, atom.z)
        b = isnothing(color_by) ? beta(atom) : color_by[i]
        if !isnothing(color_by)
            vmin, vmax = min(vmin, b), max(vmax, b)
        end
        println(cif,
            atom.flag == 0 ? "ATOM" : "HETATM", " ", natoms, " ", _cif_element(atom), " ",
            _cif_value(name(atom)), " ", _cif_value(resname(atom)), " ", _cif_value(chain(atom)), " ",
            resnum(atom), " ", @sprintf("%.3f %.3f %.3f %.2f", atom.x, atom.y, atom.z, b),
        )
        if isprotein(atom) || isnucleoside(atom)
            push!(chains, chain(atom))
            npolymer += 1
            name(atom) in ("CA", "P") && (nbackbone += 1)
            print(groups, '0')
        elseif iswater(atom)
            print(groups, '2')
        else
            nother += 1
            print(groups, '1')
        end
    end
    natoms == 0 && throw(ArgumentError("No atoms to visualize."))
    geometric_center ./= natoms
    if !isnothing(color_range)
        vmin, vmax = color_range
    end
    # 3Dmol.js gradients are not defined for an empty range
    if vmin == vmax
        vmin, vmax = vmin - 0.5, vmax + 0.5
    end
    # Color specifications
    polymer_color, other_color = if !isnothing(color_by)
        gradient = Dict("prop" => "b", "gradient" => string(colormap), "min" => vmin, "max" => vmax)
        Dict("colorscheme" => gradient), Dict("colorscheme" => gradient)
    elseif color isa AbstractString
        Dict("color" => color), Dict("color" => color)
    elseif color == :auto
        if style == :cartoon
            length(chains) > 1 ? Dict("colorscheme" => "chain") : Dict("color" => "spectrum"),
            Dict("colorscheme" => "greenCarbon")
        else
            Dict("colorscheme" => "default"), Dict("colorscheme" => "default")
        end
    else
        scheme = Dict(:chain => "chain", :ss => "ssPyMol", :element => "default", :residue => "amino")
        if color == :spectrum
            Dict("color" => "spectrum"), Dict("colorscheme" => "greenCarbon")
        elseif color == :ss
            Dict("colorscheme" => scheme[color]), Dict("colorscheme" => "greenCarbon")
        elseif haskey(scheme, color)
            Dict("colorscheme" => scheme[color]), Dict("colorscheme" => scheme[color])
        else
            throw(ArgumentError("color must be :auto, :chain, :ss, :spectrum, :element, :residue, or a color string. Got: :$color"))
        end
    end
    _style(s, color) = if s == :cartoon
        Dict("cartoon" => color)
    elseif s == :sticks
        Dict("stick" => merge(Dict("radius" => 0.2), color))
    elseif s == :ballandstick
        Dict("stick" => merge(Dict("radius" => 0.15), color), "sphere" => merge(Dict("scale" => 0.25), color))
    elseif s == :spheres
        Dict("sphere" => color)
    elseif s == :lines
        Dict("line" => color)
    end
    # Cartoons are not drawn for less than two residues: use sticks instead
    polymer_style = _style(style == :cartoon && nbackbone < 2 ? :sticks : style, polymer_color)
    other_style = _style(style == :cartoon ? :sticks : style, other_color)
    # Water is shown as lines in the presence of other molecules, otherwise with the style of other molecules
    only_water = npolymer == 0 && (nother == 0 || !ligands)
    show_water = water == :auto ? only_water : water
    water_style = only_water ? other_style : Dict("line" => Dict("colorscheme" => "default"))
    js_str(x) = replace(JSON.json(x), "</" => "<\\/")
    js = IOBuffer()
    println(js, "const cif = ", js_str(String(take!(cif))), ";")
    println(js, "const groups = ", js_str(String(take!(groups))), ";")
    println(js, "const inGroup = (g) => ({predicate: (a) => groups.charCodeAt(a.index) === 48 + g});")
    println(js, "const viewer = \$3Dmol.createViewer(document.getElementById('viewer'), {backgroundColor: ", js_str(background), "});")
    println(js, "viewer.addModel(cif, 'cif');")
    println(js, "viewer.setStyle({}, {});")
    println(js, "viewer.setStyle(inGroup(0), ", js_str(polymer_style), ");")
    ligands && println(js, "viewer.setStyle(inGroup(1), ", js_str(other_style), ");")
    show_water && println(js, "viewer.setStyle(inGroup(2), ", js_str(water_style), ");")
    # Surfaces are computed synchronously, because the 3Dmol.js web workers may not run from
    # file:// pages or from iframes embedded in notebooks or editors
    surface && println(js, "\$3Dmol.SurfaceWorker = undefined;")
    surface && println(js, "viewer.addSurface(\$3Dmol.SurfaceType.SES, {opacity: 0.6, color: 'white'}, inGroup(0));")
    if !isnothing(unitcell)
        m = unitcell isa AbstractVector ? [unitcell[1] 0 0; 0 unitcell[2] 0; 0 0 unitcell[3]] : unitcell
        a, b, c = m[:, 1], m[:, 2], m[:, 3]
        origin = unitcell_origin == :center ? geometric_center - (a + b + c) / 2 : collect(unitcell_origin)
        corner(i, j, k) = origin + i * a + j * b + k * c
        xyz(v) = Dict("x" => v[1], "y" => v[2], "z" => v[3])
        radius = max(0.05, 0.002 * maximum(norm, (a, b, c)))
        # The 12 edges connect corners that differ in one lattice vector
        for (i, j, k) in Iterators.product(0:1, 0:1, 0:1), (di, dj, dk) in ((1, 0, 0), (0, 1, 0), (0, 0, 1))
            (i + di > 1 || j + dj > 1 || k + dk > 1) && continue
            edge = Dict(
                "start" => xyz(corner(i, j, k)), "end" => xyz(corner(i + di, j + dj, k + dk)),
                "radius" => radius, "color" => unitcell_color, "fromCap" => 1, "toCap" => 1,
            )
            println(js, "viewer.addCylinder(", js_str(edge), ");")
        end
    end
    if hover
        println(js, """
        viewer.setHoverable({}, true,
            function (atom, viewer) {
                if (!atom.label) {
                    atom.label = viewer.addLabel(atom.resn + atom.resi + ":" + atom.chain + " " + atom.atom,
                        {position: atom, backgroundColor: "black", backgroundOpacity: 0.7, fontColor: "white", fontSize: 12});
                }
            },
            function (atom, viewer) {
                if (atom.label) { viewer.removeLabel(atom.label); delete atom.label; }
            });""")
    end
    println(js, "viewer.zoomTo();")
    println(js, "viewer.render();")
    html = """
    <!DOCTYPE html>
    <html>
    <head>
    <meta charset="utf-8">
    <title>PDBTools.jl - $natoms atoms</title>
    <script src="$_3DMOL_URL"></script>
    <style>html, body { margin: 0; height: 100%; overflow: hidden; } #viewer { position: relative; width: 100%; height: 100%; }</style>
    </head>
    <body>
    <div id="viewer"></div>
    <script>
    $(String(take!(js)))</script>
    </body>
    </html>
    """
    _css_size(x) = x isa Integer ? "$(x)px" : String(x)
    return StructureView(html, natoms, _css_size(width), _css_size(height))
end

# Values of the mmCIF atom_site loop: empty values are written as ".", and values
# with spaces or quotes (e.g. primed nucleotide atom names) are quoted
function _cif_value(s::AbstractString)
    s = strip(s)
    isempty(s) && return "."
    if any(isspace, s) || occursin('\'', s) || occursin('"', s) || first(s) in ('_', '#', '\$', ';', '[', ']')
        return occursin('"', s) ? "'$s'" : "\"$s\""
    end
    return s
end

# Element symbol, or the first letter of the atom name if the element is unknown
function _cif_element(atom::Atom)
    element(atom) !== nothing && return element_symbol_string(atom)
    i = findfirst(isletter, name(atom))
    return isnothing(i) ? "X" : string(name(atom)[i])
end

_escape_html_attribute(s) = replace(s, "&" => "&amp;", "\"" => "&quot;", "<" => "&lt;", ">" => "&gt;")

# Embedded in an iframe, such that the view is isolated from the page it is displayed in
function Base.show(io::IO, ::MIME"text/html", v::StructureView)
    print(io, """<iframe srcdoc="$(_escape_html_attribute(v.html))" """)
    print(io, """style="width: $(v.width); height: $(v.height); border: none;"></iframe>""")
end

Base.show(io::IO, v::StructureView) = print(io, "PDBTools.StructureView($(v.natoms) atoms)")

# VS Code plot pane
Base.show(io::IO, ::MIME"juliavscode/html", v::StructureView) = print(io, v.html)

function Base.show(io::IO, ::MIME"text/plain", v::StructureView)
    print(io, chomp("""
    PDBTools.StructureView of $(v.natoms) atoms.
        Displayed in HTML-capable environments (VS Code, Pluto, Jupyter),
        or saved to an HTML file with save("view.html", view).
    """))
end

"""
    save(filename::AbstractString, view::StructureView)

Save the interactive structure view as a standalone HTML file, that can be opened in a web browser.

"""
function save(filename::AbstractString, v::StructureView)
    open(expanduser(filename), "w") do io
        print(io, v.html)
    end
    return filename
end

@testitem "visualize" begin
    using PDBTools
    atoms = read_pdb(PDBTools.TESTPDB)
    v = visualize(atoms)
    @test v isa PDBTools.StructureView
    @test v.natoms == length(atoms)
    @test occursin("3Dmol", v.html)
    @test occursin("cartoon", v.html)
    # one group character per atom
    m = match(r"const groups = \"([012]*)\"", v.html)
    @test length(m[1]) == length(atoms)
    @test count(==('0'), m[1]) == count(isprotein, atoms)
    @test count(==('2'), m[1]) == count(iswater, atoms)
    # selections
    v = visualize(atoms, "protein and residue < 10")
    @test v.natoms == length(select(atoms, "protein and residue < 10"))
    v = visualize(atoms, at -> name(at) == "CA"; style=:spheres, color=:chain)
    @test v.natoms == count(at -> name(at) == "CA", atoms)
    @test occursin("sphere", v.html)
    # styles and colors
    for style in (:cartoon, :sticks, :ballandstick, :spheres, :lines)
        @test visualize(atoms, "protein"; style) isa PDBTools.StructureView
    end
    for color in (:auto, :chain, :ss, :spectrum, :element, :residue, "red")
        @test visualize(atoms, "protein"; color) isa PDBTools.StructureView
    end
    v = visualize(atoms; color_by=Float64.(1:length(atoms)), colormap=:roygb)
    @test occursin("roygb", v.html)
    v = visualize(atoms, "protein"; color_by=zeros(length(atoms)))
    @test occursin("\"min\":-0.5", v.html) && occursin("\"max\":0.5", v.html)
    v = visualize(atoms, "protein"; color_by=Float64.(1:length(atoms)), color_range=(0, 10))
    @test occursin("\"max\":10", v.html)
    v = visualize(atoms, "protein"; surface=true, water=true, hover=false, width=600, height=300)
    @test occursin("SurfaceType", v.html)
    @test !occursin("setHoverable", v.html)
    @test v.width == "600px" && v.height == "300px"
    @test_throws ArgumentError visualize(atoms; style=:ribbon)
    @test_throws ArgumentError visualize(atoms; color=:rainbow)
    @test_throws ArgumentError visualize(atoms; colormap=:viridis)
    @test_throws ArgumentError visualize(atoms; water=:yes)
    @test_throws ArgumentError visualize(atoms; color_by=[1.0])
    @test_throws ArgumentError visualize(atoms, "name XXXX")
    # water is shown by default only if there is nothing else to show
    @test !occursin("inGroup(2)", visualize(atoms).html)
    @test occursin("inGroup(2), {\"stick\"", visualize(atoms, "water").html)
    @test occursin("inGroup(2), {\"line\"", visualize(atoms; water=true).html)
    @test !occursin("inGroup(2)", visualize(atoms, "water"; water=false).html)
    # unit cell
    v = visualize(atoms; unitcell=[10, 20, 30])
    @test count("addCylinder", v.html) == 12
    @test occursin("\"z\":30", v.html) && !occursin("\"z\":-15", v.html)
    v = visualize(atoms; unitcell=[10, 20, 30], unitcell_origin=[0, 0, -15], unitcell_color="red")
    @test occursin("\"z\":-15", v.html) && occursin("\"red\"", v.html)
    uc = read_unitcell(PDBTools.TESTPBC)
    v = visualize(read_pdb(PDBTools.TESTPBC); unitcell=uc, unitcell_origin=:center)
    @test count("addCylinder", v.html) == 12
    @test !occursin("addCylinder", visualize(atoms).html)
    @test_throws ArgumentError visualize(atoms; unitcell=[1, 2])
    @test_throws ArgumentError visualize(atoms; unitcell=[1, 2, 3], unitcell_origin=:middle)
    @test_throws ArgumentError visualize(atoms; unitcell=[1, 2, 3], unitcell_origin=[0, 0])
    # mmCIF fields: long chain names, and quoted atom names
    atoms_long_chain = read_mmcif(PDBTools.LONG_CHAIN_STRING_CIF)
    v = visualize(atoms_long_chain)
    @test v.natoms == length(atoms_long_chain)
    @test occursin("inGroup(0), {\"stick\"", v.html) # single residue: no cartoon
    @test occursin("inGroup(0), {\"cartoon\"", visualize(atoms, "protein").html)
    long_chain = first(filter(c -> length(c) > 1, chain.(atoms_long_chain)))
    @test occursin(" $long_chain ", v.html)
    @test PDBTools._cif_value("O5'") == "\"O5'\""
    @test PDBTools._cif_value("A B") == "\"A B\""
    @test PDBTools._cif_value("") == "."
    @test PDBTools._cif_value("CA") == "CA"
    @test PDBTools._cif_element(Atom(name="CA", pdb_element="C")) == "C"
    @test PDBTools._cif_element(Atom(name="1XX", pdb_element="")) == "X"
    # display
    html = sprint(show, MIME"text/html"(), visualize(atoms, "protein"))
    @test startswith(html, "<iframe srcdoc=")
    @test !occursin("<script", html) # escaped inside srcdoc
    @test occursin("StructureView", sprint(show, MIME"text/plain"(), v))
    @test sprint(show, MIME"juliavscode/html"(), v) == v.html
    @test sprint(show, v) == "PDBTools.StructureView($(v.natoms) atoms)"
    # save
    tmpfile = tempname() * ".html"
    @test save(tmpfile, v) == tmpfile
    @test read(tmpfile, String) == v.html
end
