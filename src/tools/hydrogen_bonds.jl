using CellListMap: CellListMap, map_pairwise!, ParticleSystem1, ParticleSystem2

struct HPolarBonds
    D::Vector{Int32} # Hydrogen-bond donnor
    H::Vector{Int32} # polar hydrogen
end
CellListMap.copy_output(x::HPolarBonds) = HPolarBonds(copy(x.D), copy(x.H))
function CellListMap.reset_output!(x::HPolarBonds)
    empty!(x.D)
    empty!(x.H)
    return x
end
function CellListMap.reducer!(x::HPolarBonds, y::HPolarBonds)
    append!(x.D, y.D)
    append!(x.H, y.H)
    return x
end

@kwdef struct HBonds
    D::Vector{Int32} = Int32[] # Hydrogen bond donnor
    H::Vector{Int32} = Int32[] # polar hydrogen
    A::Vector{Int32} = Int32[] # hydrogen bond acceptor
    r::Vector{Float32} = Float32[] # distance D-A
    ang::Vector{Float32} = Float32[] # angle H-D-A
end
Base.length(x::HBonds) = length(x.D)
Base.getindex(x::HBonds, i::Integer) = (D=x.D[i], H=x.H[i], A=x.A[i], r=x.r[i], ang=x.ang[i])
function Base.show(io::IO, ::MIME"text/plain", hb::HBonds)
    print(io, chomp("""
    HBonds data structure with $(length(hb)) hydrogen-bonds.
        First hbond: (D-H---A) = $(length(hb) > 0 ? hb[1] : nothing)
        Last hbond: (D-H---A) = $(length(hb) > 0 ? hb[length(hb)] : nothing)
        - r is the distance between Donnor and Acceptor atoms (D-A)
        - ang is the angle (degrees) between H-D and A-D.
    """))
end

CellListMap.copy_output(x::HBonds) = HBonds(copy(x.D), copy(x.H), copy(x.A), copy(x.r), copy(x.ang))
function CellListMap.reset_output!(x::HBonds)
    empty!(x.D)
    empty!(x.H)
    empty!(x.A)
    empty!(x.r)
    empty!(x.ang)
    return x
end
function CellListMap.reducer!(x::HBonds, y::HBonds)
    append!(x.D, y.D)
    append!(x.H, y.H)
    append!(x.A, y.A)
    append!(x.r, y.r)
    append!(x.ang, y.ang)
    return x
end

function hbond_angle(D, H, A)
    v1 = H - D
    v2 = A - D
    return acosd(dot(v1, v2) / (norm(v1) * norm(v2)))
end

function push_hbond!(hbonds, i, j, x, y, polar_bonds, positions, unitcell, ang, ats_sel1, d2)
    # Find if i is a donnor
    ii = searchsortedfirst(polar_bonds.D, i)
    ii > length(polar_bonds.D) && return nothing
    # Might have more than one polar hydrogen
    while polar_bonds.D[ii] == i
        iH = polar_bonds.H[ii]
        xH = isnothing(unitcell) ? positions[iH] : wrap(positions[iH], x, unitcell)
        hbond_ang = hbond_angle(x, xH, y)
        if hbond_ang <= ang
            push!(hbonds.D, index(ats_sel1[i]))
            push!(hbonds.H, index(ats_sel1[iH]))
            push!(hbonds.A, index(ats_sel1[j]))
            push!(hbonds.r, sqrt(d2))
            push!(hbonds.ang, hbond_ang)
        end
        ii += 1
        ii > length(polar_bonds.D) && break
    end
    return nothing
end

function push_hbond2!(hbonds,
    i, j, x, y, polar_bonds, positions,
    unitcell, ang, ats_sel1, ats_sel2, d2
)
    # Find if i is a donnor
    ii = searchsortedfirst(polar_bonds.D, i)
    ii > length(polar_bonds.D) && return nothing
    # Might have more than one polar hydrogen
    while polar_bonds.D[ii] == i
        iH = polar_bonds.H[ii]
        xH = isnothing(unitcell) ? positions[iH] : wrap(positions[iH], x, unitcell)
        hbond_ang = hbond_angle(x, xH, y)
        if hbond_ang <= ang
            push!(hbonds.D, index(ats_sel1[i]))
            push!(hbonds.H, index(ats_sel1[iH]))
            push!(hbonds.A, index(ats_sel2[j]))
            push!(hbonds.r, sqrt(d2))
            push!(hbonds.ang, hbond_ang)
        end
        ii += 1
        ii > length(polar_bonds.D) && break
    end
    return nothing
end

function find_hbond_donnors(
    ats_sel;
    positions,
    unitcell,
    cutoff,
    parallel,
    electronegative_elements,
    d_covalent_bond,
)
    sys = ParticleSystem(;
        positions,
        unitcell,
        cutoff,
        output=HPolarBonds(Int32[], Int32[]),
        parallel,
        output_name=:polar_bonds,
    )
    polar_bonds = map_pairwise!(sys) do _, _, i, j, d2, polar_bonds
        at_i = ats_sel[i]
        at_j = ats_sel[j]
        el_i = element(at_i)
        el_j = element(at_j)
        if (el_i in electronegative_elements) & (el_j == "H")
            D = i
            H = j
        elseif (el_j in electronegative_elements) & (el_i == "H")
            D = j
            H = i
        else
            return polar_bonds
        end
        if d2 < (d_covalent_bond)^2
            push!(polar_bonds.D, D)
            push!(polar_bonds.H, H)
        end
        return polar_bonds
    end
    iord = sortperm(polar_bonds.D)
    polar_bonds.D .= polar_bonds.D[iord]
    polar_bonds.H .= polar_bonds.H[iord]
    return polar_bonds
end

#=

    process_selections(selections) -> Vector{Pair{String,String}}

Convert input selections into a list of pairs. If no selection is provided, returns ["all"=>"all"].

=#
function process_selections(selections)
    if isnothing(first(selections))
        return ["all" => "all"]
    end
    return [sel isa String ? sel => sel : first(sel) => last(sel) for sel in selections]
end

struct SelectionData{A<:PDBTools.Atom,}
    ats::Vector{A}
    inds::Vector{Int32}
    polar_bonds::HPolarBonds
end
_key_name(sel1, sel2) = "$sel1 => $sel2"


#=

    initialize_hbonds_data(sim, selection_pairs)

Initialize the results dictionary and process selection data for each unique selection.

=#
function initialize_hbonds_data(
    atoms::AbstractVector{<:PDBTools.Atom},
    selection_pairs::Vector{Pair{String,String}};
    unitcell=nothing,
    parallel::Bool=true,
    donnor_acceptor_distance::Real=3.5f0,
    electronegative_elements=("N", "O", "F", "S"),
    d_covalent_bond::Real=1.2f0,
)
    selection_data = Dict{String,SelectionData}()
    for selection_pair in selection_pairs
        for sel in selection_pair
            if !haskey(selection_data, sel)
                ats_sel = select(atoms, sel)
                inds_sel = index.(ats_sel)
                polar_bonds = find_hbond_donnors(
                    ats_sel,
                    positions=coor(ats_sel),
                    unitcell=unitcell,
                    cutoff=donnor_acceptor_distance,
                    parallel=parallel,
                    electronegative_elements=electronegative_elements,
                    d_covalent_bond=d_covalent_bond,
                )
                selection_data[sel] = SelectionData(ats_sel, inds_sel, polar_bonds)
            end
        end
    end
    return selection_data
end

#=

Initialize CellListMap particle systems for each selection pair.

=#
function setup_particle_systems(
    selection_pairs,
    selection_data,
    unitcell,
    donnor_acceptor_distance,
    parallel
)
    systems = Dict{String,Union{CellListMap.ParticleSystem1,CellListMap.ParticleSystem2}}()
    for selection_pair in selection_pairs
        sel1, sel2 = first(selection_pair), last(selection_pair)
        key = _key_name(sel1, sel2)
        if sel1 == sel2
            s1 = selection_data[sel1]
            systems[key] = ParticleSystem(
                xpositions=coor.(s1.ats),
                unitcell=unitcell,
                cutoff=donnor_acceptor_distance,
                output=HBonds(),
                parallel=parallel,
                output_name=:hydrogen_bonds
            )
        else
            s1 = selection_data[sel1]
            s2 = selection_data[sel2]
            systems[key] = ParticleSystem(
                xpositions=coor.(s1.ats),
                ypositions=coor.(s2.ats),
                unitcell=unitcell,
                cutoff=donnor_acceptor_distance,
                output=HBonds(),
                parallel=parallel,
                output_name=:hydrogen_bonds
            )
        end
    end
    return systems
end

function compute_hbonds!(sys::ParticleSystem1, s1, angle_cutoff, electronegative_elements)
    map_pairwise!(sys) do x, y, i, j, d2, hbonds
        el_i = element(s1_ats[i])
        el_j = element(s2.ats[j])
        if (el_i in electronegative_elements) & (el_j in electronegative_elements)
            push_hbond!(hbonds,
                i, j, x, y, s1.polar_bonds, sys.xpositions,
                sys.unitcell, angle_cutoff, s1.ats, d2
            )
            push_hbond!(hbonds,
                j, i, y, x, s1.polar_bonds, sys.xpositions,
                sys.unitcell, angle_cutoff, s1.ats, d2
            )
        end
        return hbonds
    end
end

function compute_hbonds!(sys::ParticleSystem2, s1, s2, sel1, sel2, angle_cutoff, electronegative_elements)
    map_pairwise!(sys) do x, y, i, j, d2, hbonds
        at_i = s1.ats[i]
        at_j = s2.ats[j]
        if index(at_i) == index(at_j)
            throw(ArgumentError("""\n
                Different selections cannot overlap. Detected atom $(index(at_i)) in both selections \"$sel1\" and \"$sel2\".
                """))
        end
        el_i = element(at_i)
        el_j = element(at_j)
        if (el_i in electronegative_elements) & (el_j in electronegative_elements)
            push_hbond2!(hbonds,
                i, j, x, y, s1.polar_bonds, sys.xpositions,
                sys.unitcell, angle_cutoff, s1.ats, s2.ats, d2
            )
            push_hbond2!(hbonds,
                j, i, y, x, s2.polar_bonds, sys.ypositions,
                sys.unitcell, angle_cutoff, s2.ats, s1.ats, d2
            )
        end
        return hbonds
    end
end

"""
    hydrogen_bonds(atoms, sel, sel1 => sel2, ... ; kargs...)

Function to find hydrogen bonds in a set of atoms, or among two sets of atoms. The structure must 
contain Hydrogen atoms.

### Arguments

- `atoms`: Vector of atoms, or structure component (model, chain, segment, residue) to be analyzed.
and, optionally, the selections or selection pairs for which the hydrogen bonds must be computed:
- `sel::String`: Selection string, e. g. `"protein"`.
- `sel1 => sel2::Pair{String,String}`: Pair of selection strings, e. g. `"resname ARG" => "resname GLU"`.

The two selections of each pair, if different, must not have overalpping atoms (and error will the thrown).
If no selection is provided, the hydrogen bonds of the complete structure will be computed.

### Optional keyword arguments

- `unitcell::Union{Nothing,AbstractVecOrMat}=nothing`: Unit cell for periodic boundary conditions.
- `donnor_acceptor_distance::Real=3.5f0`: Maximum distance between donnor and acceptor to consider a hydrogen bond.
- `angle_cutoff::Real=30`: Maximum angle (in degrees) between donnor-hydrogen-acceptor to consider a hydrogen bond.
- `electronegative_elements=("N", "O", "F", "S")`: Elements considered electronegative for hydrogen bonding.
- `d_covalent_bond::Real=1.2f0`: Maximum distance between donnor and hydrogen to consider a covalent bond.
- `parallel::Bool=false`: Whether to use parallel computation.

### Returns

- `HBonds`: A data structure containing the found hydrogen bonds, where each element is 
  a named tuple `(D, H, A, r, ang)` where the fields correspond to the donnor, hydrogen
  and acceptor atoms, the distance between donnor and acceptor atoms, and the angle. 

# Example

```jldoctest; filter = r"(\\d*)\\.(\\d{2})\\d+" => s"\\1.\\2***"
julia> using PDBTools

julia> pdb = read_pdb(PDBTools.test_dir*"/hbonds.pdb", "model 1");

julia> uc = read_unitcell(PDBTools.test_dir*"/hbonds.pdb");

julia> hbs = hydrogen_bonds(pdb, "protein"; unitcell=uc) # Single set of atoms: selection is optional
OrderedCollections.OrderedDict{String, PDBTools.HBonds} with 1 entry:
  "protein => protein" => HBonds(Int32[1, 1, 271, 37, 1020, 237, 56, 76, 1060, 204  …  748, 813, 828, 871, 863, 877, 96…

julia> hbs["protein => protein"] # Summary
HBonds data structure with 63 hydrogen-bonds.
    First hbond: (D-H---A) = (D = 1, H = 2, A = 286, r = 2.6871147f0, ang = 10.643958f0)
    Last hbond: (D-H---A) = (D = 1014, H = 1017, A = 1032, r = 2.5816715f0, ang = 12.714139f0)
    - r is the distance between Donnor and Acceptor atoms (D-A)
    - ang is the angle (degrees) between H-D and A-D.

julia> hbs["protein => protein"][1] # first h-bond
(D = 1, H = 2, A = 286, r = 2.6871147f0, ang = 10.643958f0)

julia> hbs = hydrogen_bonds(pdb, "protein", "protein" => "resname SOL"; unitcell=uc) # Multiple selections
OrderedCollections.OrderedDict{String, PDBTools.HBonds} with 2 entries:
  "protein => protein"     => HBonds(Int32[1, 1, 271, 37, 1020, 237, 56, 76, 1060, 204  …  748, 813, 828, 871, 863, 877…
  "protein => resname SOL" => HBonds(Int32[1406, 1583, 1799, 2027, 789, 3503, 1169, 184, 3914, 4304  …  1224, 38768, 12…

julia> hbs["protein => protein"]
HBonds data structure with 63 hydrogen-bonds.
    First hbond: (D-H---A) = (D = 1, H = 2, A = 286, r = 2.6871147f0, ang = 10.643958f0)
    Last hbond: (D-H---A) = (D = 1014, H = 1017, A = 1032, r = 2.5816715f0, ang = 12.714139f0)
    - r is the distance between Donnor and Acceptor atoms (D-A)
    - ang is the angle (degrees) between H-D and A-D.

julia> hbs["protein => resname SOL"]
HBonds data structure with 138 hydrogen-bonds.
    First hbond: (D-H---A) = (D = 1406, H = 1407, A = 160, r = 2.9361732f0, ang = 6.771988f0)
    Last hbond: (D-H---A) = (D = 41798, H = 41800, A = 395, r = 2.6894214f0, ang = 10.623453f0)
    - r is the distance between Donnor and Acceptor atoms (D-A)
    - ang is the angle (degrees) between H-D and A-D.
```

!!! note
    This function does not use topology information. It identified polar hydrogens based on distance criteria only,
    where `d_covalent_bond` is the criterium for identifying covalent bonds between donnor and hydrogen atoms.

"""
function hydrogen_bonds(
    structure::Union{AbstractVector{<:Atom},AbstractStructuralElement},
    selections::Union{Nothing,String,Pair{String,String}}...=nothing;
    donnor_acceptor_distance::Real=3.5f0,
    angle_cutoff::Real=30,
    electronegative_elements=("N", "O", "F", "S"),
    unitcell::Union{Nothing,AbstractVecOrMat}=nothing,
    d_covalent_bond::Real=1.2f0,
    parallel::Bool=false,
)
    ats = get_atoms(structure)
    selection_pairs = process_selections(selections)

    # Initialize results and process selections
    selection_data = initialize_hbonds_data(
        ats, selection_pairs;
        parallel=parallel,
        donnor_acceptor_distance=donnor_acceptor_distance,
        electronegative_elements=electronegative_elements,
        d_covalent_bond=d_covalent_bond
    )

    systems = setup_particle_systems(
        selection_pairs, selection_data,
        unitcell, donnor_acceptor_distance, parallel
    )

    hbonds = OrderedDict{String,HBonds}()
    for selection_pair in selection_pairs
        sel1, sel2 = first(selection_pair), last(selection_pair)
        key = _key_name(sel1, sel2)
        sys = systems[key]
        if sel1 == sel2
            s1 = selection_data[sel1]
            compute_hbonds!(sys, s1, angle_cutoff, electronegative_elements)
        else
            s1 = selection_data[sel1]
            s2 = selection_data[sel2]
            compute_hbonds!(sys, s1, s2, sel1, sel2, angle_cutoff, electronegative_elements)
        end
        hbonds[key] = sys.hydrogen_bonds
    end
    return hbonds
end

@testitem "Hydrogen Bonds" begin
    using PDBTools
    using ShowMethodTesting
    pdb = read_pdb(PDBTools.test_dir * "/hbonds.pdb")
    models = collect(eachmodel(pdb))

    #
    # Hbonds within same set
    #
    nhb = [63, 62, 56, 57, 62] # checked with gmx hbond
    for (i, m) in enumerate(models)
        hbs = hydrogen_bonds(m, "protein")
        @test length(hbs["protein => protein"]) == nhb[i]
        hbs = hydrogen_bonds(m, "protein"; parallel=true)
        @test length(hbs["protein => protein"]) == nhb[i]
    end
    uc = read_unitcell(PDBTools.test_dir * "/hbonds.pdb")
    hbs = hydrogen_bonds(models[1], "protein"; unitcell=uc)
    @test length(hbs["protein => protein"]) == 63

    prot = select(models[1], "protein")
    hbs = hydrogen_bonds(prot)
    @test length(hbs["all => all"]) == 63

    pbcs = [
        [88.749199, 90.714630, 94.511879],
        [88.838280, 90.805679, 94.606750],
        [88.764839, 90.730606, 94.528542],
        [88.883881, 90.852287, 94.655304],
        [88.828445, 90.795631, 94.596283],
    ]
    # the following tests are using isapprox because of changes in LinearAlgebra can make them
    # fail in different versions because of the rounding associated to PBCs. In Julia 1.12 they
    # provide exact results. 
    nhb = [17964, 17915, 17945, 17977, 17852] # checked with gmx hbond
    for (i, m) in enumerate(models)
        local uc = pbcs[i]
        hbs = hydrogen_bonds(m, "resname HOH SOL"; unitcell=uc)
        @test length(hbs["resname HOH SOL => resname HOH SOL"]) ≈ nhb[i] atol = 10
        hbs = hydrogen_bonds(m, "resname HOH SOL"; unitcell=uc, parallel=true)
        @test length(hbs["resname HOH SOL => resname HOH SOL"]) ≈ nhb[i] atol = 10
    end

    #
    # Hbonds among two sets
    #
    nhb = [140, 130, 142, 135, 143] # checked with gmx hbond
    for (i, m) in enumerate(models)
        local uc = pbcs[i]
        hbs = hydrogen_bonds(m, "protein" => "resname HOH SOL")
        @test length(hbs["protein => resname HOH SOL"]) == nhb[i]
        hbs = hydrogen_bonds(m, "protein" => "resname HOH SOL"; unitcell=uc)
        @test length(hbs["protein => resname HOH SOL"]) == nhb[i]
        hbs = hydrogen_bonds(m, "protein" => "resname HOH SOL"; unitcell=uc, parallel=true)
        @test length(hbs["protein => resname HOH SOL"]) == nhb[i]
    end
    nhb = [156, 143, 150, 141, 156] # checked with gmx hbond
    for (i, m) in enumerate(models)
        local uc = pbcs[i]
        hbs = hydrogen_bonds(m, "resname HOH" => "resname SOL"; unitcell=uc)
        @test length(hbs["resname HOH => resname SOL"]) == nhb[i]
    end

    hbs = hydrogen_bonds(models[1], "protein"; unitcell=pbcs[1])
    @test parse_show(hbs["protein => protein"]) ≈ """
        HBonds data structure with 63 hydrogen-bonds.
            First hbond: (D-H---A) = (D = 1, H = 4, A = 267, r = 2.6454127f0, ang = 4.0603805f0)
            Last hbond: (D-H---A) = (D = 1169, H = 1170, A = 619, r = 2.6004055f0, ang = 9.524022f0)
            - r is the distance between Donnor and Acceptor atoms (D-A)
            - ang is the angle (degrees) between H-D and A-D.
        """

    # Test set overlap error
    @test_throws ArgumentError hydrogen_bonds(models[1], "resname ARG" => "protein"; unitcell=pbcs[1])

end