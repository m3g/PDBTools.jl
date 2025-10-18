
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
    append!(x.D,y.D)
    append!(x.H,y.H)
    return x
end

struct HBonds
    D::Vector{Int32} # Hydrogen bond donnor
    H::Vector{Int32} # polar hydrogen
    A::Vector{Int32} # hydrogen bond acceptor
    r::Vector{Float32} # distance D-A
    ang::Vector{Float32} # angle H-D-A
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
    append!(x.D,y.D)
    append!(x.H,y.H)
    append!(x.A,y.A)
    append!(x.r,y.r)
    append!(x.ang,y.ang)
    return x
end

function hbond_angle(D,H,A)
    v1 = H - D
    v2 = A - D 
    return acosd(dot(v1,v2) / (norm(v1)*norm(v2)))
end

function push_hbond!(i, j, x, y, polar_bonds, sys, ang, hbonds, ats_sel1, d2)
    # Find if i is a donnor
    ii = searchsortedfirst(polar_bonds.D, i)
    ii > length(polar_bonds.D) && return nothing
    # Might have more than one polar hydrogen
    while polar_bonds.D[ii] == i
        iH = polar_bonds.H[ii]
        xH = wrap(sys.positions[iH], x, sys.unitcell)
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

function find_hbond_donnors(
    ats_sel;
    positions, 
    unitcell, 
    cutoff,
    parallel,
    electronegative_elements,
    d_covalent_bond,
)
    sys_polar_bonds = ParticleSystem(;
        positions,
        unitcell,
        cutoff,
        output=HPolarBonds(Int32[], Int32[]),
        parallel,
        output_name=:polar_bonds,
    )
    polar_bonds = map_pairwise!(
        (x,y,i,j,d2,polar_bonds) -> begin
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
        end,
        sys_polar_bonds,
    )
    iord = sortperm(polar_bonds.D)
    polar_bonds.D .= polar_bonds.D[iord]
    polar_bonds.H .= polar_bonds.H[iord]
    return polar_bonds
end

"""
    hydrogen_bonds(ats::AbstractVector{<:PDBTools.Atom}, sel::Union{Function, Strign}=at -> true)

!!! warning
    Experimental feature: interface changes can occur in non-breaking releases.

Function to find hydrogen bonds in a set of atoms.

### Arguments

- `ats::AbstractVector{<:PDBTools.Atom}`: Vector of atoms to analyze.
- `sel::Union{Function, String}=at -> true`: Selection of atoms to consider. Can be a function or a selection string.

### Keyword Arguments

- `donnor_acceptor_distance::Real=3.5f0`: Maximum distance between donnor and acceptor to consider a hydrogen bond.
- `angle_cutoff::Real=30`: Maximum angle (in degrees) between donnor-hydrogen-acceptor to consider a hydrogen bond.
- `electronegative_elements=("N", "O", "F", "S")`: Elements considered electronegative for hydrogen bonding.
- `unitcell::Union{Nothing,AbstractVecOrMat}=nothing`: Unit cell for periodic boundary conditions.
- `d_covalent_bond::Real=1.2f0`: Maximum distance between donnor and hydrogen to consider a covalent bond.
- `parallel::Bool=false`: Whether to use parallel computation.

### Returns

- `HBonds`: A data structure containing the found hydrogen bonds.

!!! note
    This function does not use topology information. It identified polar hydrogens based on distance criteria only,
    where `d_covalent_bond` is the criterium for identifying covalent bonds between donnor and hydrogen atoms.

"""
function hydrogen_bonds(
    ats::AbstractVector{<:PDBTools.Atom},
    sel::Union{Function, String}=at -> true;
    donnor_acceptor_distance::Real=3.5f0,
    angle_cutoff::Real=30,
    electronegative_elements=("N", "O", "F", "S"),
    unitcell::Union{Nothing,AbstractVecOrMat}=nothing,
    d_covalent_bond::Real=1.2f0,
    parallel::Bool=false,
)
    # Select desired atoms
    ats_sel = select(ats, sel)
    ats_sel_positions = coor.(ats_sel)

    #
    # Find all possible hydrogen-bonding donors
    #
    polar_bonds = find_hbond_donnors(
        ats_sel;
        positions=ats_sel_positions, 
        unitcell, 
        cutoff=donnor_acceptor_distance,
        parallel,
        electronegative_elements,
        d_covalent_bond,
    )

    #
    # Find all hydrogen bonds
    #
    sys_hbonds = ParticleSystem(
        positions=ats_sel_positions,
        unitcell=unitcell,
        cutoff=donnor_acceptor_distance,
        output=HBonds(Int32[], Int32[], Int32[], Float32[], Float32[]),
        parallel=parallel,
        output_name=:hbonds
    )
    map_pairwise!(
        (x,y,i,j,d2,hbonds) -> begin
            el_i = element(ats_sel[i])
            el_j = element(ats_sel[j])
            if (el_i in electronegative_elements) & (el_j in electronegative_elements) 
                push_hbond!(i, j, x, y, polar_bonds, sys_hbonds, angle_cutoff, hbonds, ats_sel, d2)
                push_hbond!(j, i, y, x, polar_bonds, sys_hbonds, angle_cutoff, hbonds, ats_sel, d2)
            end
            return hbonds
        end,
        sys_hbonds,
    )
    return sys_hbonds.hbonds
end

function hydrogen_bonds2(
    ats::AbstractVector{<:PDBTools.Atom},
    sel1::Union{Function, String}=at -> true,
    sel2::Union{Function, String}=at -> true;
    donnor_acceptor_distance::Real=3.5f0,
    angle_cutoff::Real=30,
    electronegative_elements=("N", "O", "F", "S"),
    unitcell::Union{Nothing,AbstractVecOrMat}=nothing,
    d_covalent_bond::Real=1.2f0,
    parallel::Bool=false,
)
    # Select desired atoms
    ats_sel1 = select(ats, sel1)
    ats_sel2 = select(ats, sel2)
    ats_sel1_positions = coor.(ats_sel1)
    ats_sel2_positions = coor.(ats_sel2)

    #
    # Find all possible hydrogen-bonding donors
    #
    polar_bonds1 = find_hbond_donnors(
        ats_sel1;
        positions=ats_sel1_positions, 
        unitcell, 
        cutoff=donnor_acceptor_distance,
        parallel,
        electronegative_elements,
        d_covalent_bond,
    )
    polar_bonds2 = find_hbond_donnors(
        ats_sel2;
        positions=ats_sel2_positions,
        unitcell,
        cutoff=donnor_acceptor_distance,
        parallel,
        electronegative_elements,
        d_covalent_bond,
    )


    #
    # Now, find all hydrogen bonds
    #
    sys_hbonds = ParticleSystem(
        positions=ats_sel1_positions,
        unitcell=unitcell,
        cutoff=donnor_acceptor_distance,
        output=HBonds(Int32[], Int32[], Int32[], Float32[], Float32[]),
        parallel=parallel,
        output_name=:hbonds
    )
    map_pairwise!(
        (x,y,i,j,d2,hbonds) -> begin
            el_i = element(ats_sel1[i])
            el_j = element(ats_sel1[j])
            if (el_i in electronegative_elements) & (el_j in electronegative_elements) 
                push_hbond!(i, j, x, y, polar_bonds, sys_hbonds, angle_cutoff, hbonds, ats_sel1, d2)
                push_hbond!(j, i, y, x, polar_bonds, sys_hbonds, angle_cutoff, hbonds, ats_sel1, d2)
            end
            return hbonds
        end,
        sys_hbonds,
    )
    return sys_hbonds.hbonds
end

@testitem "Hydrogen Bonds" begin
    using PDBTools
    pdb = read_pdb(PDBTools.test_dir*"/hbonds.pdb")
    models = collect(eachmodel(pdb))
    #
    # Single structure
    #
    nhb = [63, 62, 56, 57, 62] # checked with gmx hbond
    for (i, m) in enumerate(models)
        @test length(hydrogen_bonds(m, "protein")) == nhb[i]
        @test length(hydrogen_bonds(m, "protein"; parallel=true)) == nhb[i]
    end
    uc = read_unitcell(PDBTools.test_dir*"/hbonds.pdb")
    @test length(hydrogen_bonds(models[1], "protein"; unitcell=uc)) == 63

    pbcs = [
        [ 88.749199, 90.714630, 94.511879],
        [ 88.838280, 90.805679, 94.606750],
        [ 88.764839, 90.730606, 94.528542],
        [ 88.883881, 90.852287, 94.655304],
        [ 88.828445, 90.795631, 94.596283],
    ]
    nhb = [17964, 17915, 17945, 17977, 17852] # checked with gmx hbond
    for (i, m) in enumerate(models)
        uc = pbcs[i]
        @test length(hydrogen_bonds(m, "resname HOH SOL"; unitcell=uc)) == nhb[i]
        @test length(hydrogen_bonds(m, "resname HOH SOL"; unitcell=uc, parallel=true)) == nhb[i]
    end

end