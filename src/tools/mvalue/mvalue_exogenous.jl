#=

The following functions are used to compute the m-values using less integrated SASA calculations,
which can be used for testing purposes, or for direct comparison with published results. 

These functions are considered public but are not exported.

=#
"""
    mvalue_delta_sasa(; model=MoeserHorinek, cosolvent="urea", atoms:AbstractVector{<:PDBTools.Atom}, sasas, type=1)

Calculates the m-value (transfer free energy of a protein in 1M solution) using the Tanford transfer model,
as implemented by Moeser and Horinek [1] or by Auton and Bolen [2,3].

# Arguments

- `model`: The model to be used. Must be `MoeserHorinek` or `AutonBolen`. `MoeserHorinek` is only implemented for `cosolvent="urea"`,
   and should be more precise in that case. Other solvents are available for `AutonBolen`.
- `cosolvent::String`: One of $(join('"' .* sort!(unique(keys(PDBTools.cosolvent_column)) .* '"'; by=lowercase),", "))
- `atoms::AbstractVector{<:PDBTools.Atom}`: Vector containing the atoms of the structure.
- `sasas::AbstractDict{String, AbstractDict{Symbol, Float64}}`: A dictionary containing the change in solvent accessible surface area (SASA)
  upon denaturation for each amino acid type. This data can be obtained from the `creamer_delta_sasa` function, the m-value server, or calculated using GROMACS:
    - The `creamer_delta_sasa` function provides estimated variations in SASA of a protein, using the Creamer unfolded model.
    - The output of the server can be parsed using the `parse_mvalue_server_sasa` function defined in this module.
    - Compute the SASA with `delta_sasa_per_restype`, a SASA calculation utility implemented in PDBTools.jl.
    - SASA values can be calculated using GROMACS with the `gmx_delta_sasa_per_restype` function defined in this module.

- `type::Int`: Specifies which SASA value to use from the provided data, because the server provides minimum, average, and maximum values,
    according to different denatured models for the protein. The recommended value is `2` for comparison with experimental data.
    Normally, GROMACS calculations will provide a single value, so `type=1` should be used in that case.

# Returns

A named tuple with the following fields:
- `tot`: Total transfer free energy (kcal/mol).
- `bb`: Contribution from the backbone (kcal/mol).
- `sc`: Contribution from the side chains (kcal/mol).
- `restype`: A dictionary with the transfer free energy contributions per residue type (kcal/mol).

Each entry in the dictionary is a named tuple with `bb` and `sc` fields representing the backbone and side chain contributions, respectively.

# Example calls

```julia
using PDBTools
using PDBTools: mvalue_delta_sasa,
                delta_sasa_per_restype,
                creamer_delta_sasa,
                parse_mvalue_server_sasa,
                gmx_delta_sasa_per_restype,

protein = read_pdb("protein.pdb")

# Using SASA values calculated with PDBTools.jl
sasas=delta_sasa_per_restype(native=read_pdb("native.pdb"), desnat=read_pdb("desnat.pdb"))
mvalue_delta_sasa(; model=AutonBolen, cosolvent="TMAO", atoms=protein, sasas=sasas)

# Using SASA values computed for the Creamer denatured states
sasas_from_creamer=creamer_delta_sasa(protein)
mvalue_delta_sasa(; model=MoeserHorinek, cosolvent="urea", atoms=protein, sasas=sasas_from_creamer, type=2)

# Using SASA values from the m-value server
sasas_from_server=parse_mvalue_server_sasa(server_output)
mvalue_delta_sasa(; model=MoeserHorinek, cosolvent="urea", atoms=protein, sasas=sasas_from_server, type=2)

# Using SASA values calculated with GROMACS
sasas_gmx=gmx_delta_sasa_per_restype(native_pdb="native.pdb", desnat_pdb="desnat.pdb")
mvalue_delta_sasa(; model=AutonBolen, cosolvent="TMAO", atoms=protein, sasas=sasas_gmx)
```

## References

1. https://doi.org/10.1021/acs.jpcb.7b02138
2. https://doi.org/10.1016/s0076-6879(07)28023-1
3. https://www.pnas.org/doi/10.1073/pnas.0706251104

"""
function mvalue_delta_sasa(;
    model::Type{<:MValueModel}=MoeserHorinek,
    cosolvent::String="urea",
    atoms::AbstractVector{<:PDBTools.Atom}, sasas, type=1
)
    protein = select(atoms, "protein")
    residue_names = unique(resname.(protein))
    DeltaG_per_residue = OrderedDict{String,@NamedTuple{bb::Float64, sc::Float64}}()
    for rname in residue_names
        rtype = threeletter(rname) # convert non-standard residue names in types (e. g. HSD -> HIS)
        bb_type, sc_type = tfe_asa(model, cosolvent, rtype)
        DeltaG_per_residue[rtype] = (
            bb=bb_type * sasas[rtype][:bb][type],
            sc=sc_type * sasas[rtype][:sc][type]
        )
    end
    DeltaG_BB = sum(getfield(DeltaG_per_residue[key], :bb) for key in keys(DeltaG_per_residue))
    DeltaG_SC = sum(getfield(DeltaG_per_residue[key], :sc) for key in keys(DeltaG_per_residue))
    DeltaG = DeltaG_BB + DeltaG_SC
    return (tot=DeltaG, bb=DeltaG_BB, sc=DeltaG_SC, restype=DeltaG_per_residue)
end

"""
    delta_sasa_per_restype(; 
        native::AbstractVector{<:PDBTools.Atom}, 
        desnat::AbstractVector{<:PDBTools.Atom}
    )

Calculates the change in solvent accessible surface area (SASA) upon denaturation for each amino acid type
using PDBTools. Returns a dictionary that can be directly used as input to the `mvalue` function.

# Arguments

- `native`: Vector of PDBTools.Atom objects for the native structure.
- `desnat`: Vector of PDBTools.Atom objects for the denatured structure.

# Returns

A dictionary where each key is an amino acid three-letter code (e.g., "ALA", "PHE"), and the value
is another dictionary with two keys: `:sc` for side chain SASA values and `:bb` for backbone SASA values.
Each of these keys maps to a tuple containing a single Float64 value representing the change in SASA upon denaturation in Å².

# Optional arguments

- `n_dots::Int=500`: Sets the precision of the SASA calculation (greater is better).
- `backbone::Function = at -> name(at) in ("N", "CA", "C", "O")`: Define what is a backbone atom.
- `sidechain::Function = at -> !(name(at) in ("N", "CA", "C", "O"))`: Define what is a sidechain atom.
- `ignore_hydrogen::Bool = true`: By default, ignore all Hydrogen atoms of the structure.
- `unitcell=nothing`: By default, do not use periodic boundary conditions. To use PBCs, define A
  unitcell by providing either a 3x3 matrix or, for orthorhombic cells, a vector of length 3 of cell sides.

"""
function delta_sasa_per_restype(;
    native::AbstractVector{<:PDBTools.Atom},
    desnat::AbstractVector{<:PDBTools.Atom},
    n_dots::Int=500,
    backbone::Function=at -> name(at) in ("N", "CA", "C", "O"),
    sidechain::Function=at -> !(name(at) in ("N", "CA", "C", "O")),
    ignore_hydrogen::Bool=true,
    unitcell=nothing,
)
    keepH(at) = ignore_hydrogen ? !(element(at) == "H") : true
    native_atoms = select(native, at -> isprotein(at) & keepH(at))
    desnat_atoms = select(desnat, at -> isprotein(at) & keepH(at))
    native_atomic_sasa = PDBTools.sasa_particles(native_atoms; n_dots, unitcell)
    desnat_atomic_sasa = PDBTools.sasa_particles(desnat_atoms; n_dots, unitcell)
    sasas = OrderedDict{String,OrderedDict{Symbol,Float32}}()
    for rname in unique(resname.(eachresidue(native_atoms)))
        bb_native = PDBTools.sasa(native_atomic_sasa, at -> resname(at) == rname && backbone(at))
        bb_desnat = PDBTools.sasa(desnat_atomic_sasa, at -> resname(at) == rname && backbone(at))
        if !ignore_hydrogen || rname != "GLY"
            sc_native = PDBTools.sasa(native_atomic_sasa, at -> resname(at) == rname && sidechain(at))
            sc_desnat = PDBTools.sasa(desnat_atomic_sasa, at -> resname(at) == rname && sidechain(at))
        else
            sc_native = 0.0
            sc_desnat = 0.0
        end
        # convert to standard residue (e. g. HSD -> HIS) name and save
        sasas[threeletter(rname)] = OrderedDict(:sc => sc_desnat - sc_native, :bb => bb_desnat - bb_native)
    end
    return sasas # Å^2  
end

"""
    parse_mvalue_server_sasa(string::AbstractString)

Parses the SASA output from the m-value calculator server (http://best.bio.jhu.edu/mvalue/), into a dictionary
that can be directly used as input to the `mvalue` function.

The input string should contain lines formatted as follows, and correspond to the SASA values for each amino acid type:

```julia
sasa_from_server = \"\"\"
ALA 	8 	 (    11.1)     79.1 [   147.1] | (   -13.0)     51.4 [   115.8] 
PHE 	3 	 (   166.9)    197.1 [   230.2] | (    29.4)     56.4 [    83.4] 
LEU 	7 	 (   475.2)    532.2 [   589.3] | (    89.3)    145.3 [   201.3] 
...
LYS 	6 	 (   171.5)    220.4 [   269.3] | (    -4.5)     42.0 [    88.5] 
ARG 	1 	 (   110.2)    124.4 [   138.6] | (    17.1)     25.0 [    33.0] 
CYS 	0 	 (     0.0)      0.0 [     0.0] | (     0.0)      0.0 [     0.0] 
\"\"\"
```

This data can be found in the output of the server, under the title 
`"Sidechain and Backbone changes in Accessible Surface Area"`.

The function returns a dictionary where each key is an amino acid three-letter code (e.g., `"ALA"`, `"PHE"`), and the value 
is another dictionary with two keys: `:sc` for side chain SASA values and `:bb` for backbone SASA values. 
Each of these keys maps to a tuple containing three Float64 values representing the minimum, average, and maximum SASA values in Å².

"""
function parse_mvalue_server_sasa(string::AbstractString)
    sasa = OrderedDict{String,OrderedDict{Symbol,Tuple{Float64,Float64,Float64}}}()
    for line in split(string, "\n")
        # Replace (, ) and [, ], | with spaces
        line = replace(line, '(' => ' ', ')' => ' ', '[' => ' ', ']' => ' ', '|' => ' ')
        data = split(line)
        if length(data) < 8
            continue
        end
        resname = data[1]
        sc1 = parse(Float64, data[3])
        sc2 = parse(Float64, data[4])
        sc3 = parse(Float64, data[5])
        bb1 = parse(Float64, data[6])
        bb2 = parse(Float64, data[7])
        bb3 = parse(Float64, data[8])
        sasa[resname] = OrderedDict(:sc => (sc1, sc2, sc3), :bb => (bb1, bb2, bb3))
    end
    return sasa
end

#=
    read_gmx_delta_sasa_per_restype_values(filename::String, n)

Reads the output of `gmx sasa` and returns the SASA values.
`n` is the number of surfaces calculated (1 for BB only, 2 for SC and BB, for example).

=#
function read_gmx_delta_sasa_per_restype_values(filename::String, n)
    local sasa_values
    open(filename, "r") do io
        for line in eachline(io)
            if startswith(line, r"@|#")
                continue
            end
            sasa_values = parse.(Float64, split(line)[3:2+n])
        end
    end
    return sasa_values
end

"""
    gmx_delta_sasa_per_restype(; native_pdb::AbstractString, desnat_pdb::AbstractString)

Calculates the change in solvent accessible surface area (SASA) upon denaturation for each amino acid type
using GROMACS. Returns a dictionary that can be directly used as input to the `mvalue` function.

!!! note
    This function requires GROMACS (`gmx sasa` executable) to be installed and accessible from the command line.
    The path to the `gmx` executable can be provided with the `gmx` keyword.

# Arguments

- `native_pdb::AbstractString`: Path to the PDB file of the native protein structure.
- `desnat_pdb::AbstractString`: Path to the PDB file of the denatured protein structure.

# Optional arguments 

- `gmx`: the path to the `gmx` GROMACS exectuable (by default it expects `gmx` to be on the path).
- `n_dots::Int`: sets the precision of the SASA grid (greater is better).
- `backbone::Function = at -> name(at) in ("N", "CA", "C", "O")`: Define what is a backbone atom.
- `sidechain::Function = at -> !(name(at) in ("N", "CA", "C", "O"))`: Define what is a sidechain atom.
- `ignore_hydrogen::Bool=true`: By default, ignore all hydrogen atoms.

# Returns

A dictionary where each key is an amino acid three-letter code (e.g., "ALA", "PHE"), and the value
is another dictionary with two keys: `:sc` for side chain SASA values and `:bb` for backbone SASA values.
Each of these keys maps to a tuple containing a single Float64 value representing the change in SASA upon denaturation in Å².

"""
function gmx_delta_sasa_per_restype(;
    native_pdb::AbstractString,
    desnat_pdb::AbstractString,
    n_dots::Int=500,
    ignore_hydrogen::Bool=true,
    gmx="gmx",
    backbone::Function=at -> name(at) in ("N", "CA", "C", "O"),
    sidechain::Function=at -> !(name(at) in ("N", "CA", "C", "O")),
)
    sasas = OrderedDict{String,OrderedDict{Symbol,Float64}}()
    p = read_pdb(native_pdb, "protein")
    for rname in unique(resname.(eachresidue(p)))
        sasa_bb_native, sasa_sc_native = 100 .* gmx_sasa_per_restype(native_pdb, rname; gmx, n_dots, backbone, sidechain, ignore_hydrogen) # returns in nm^2
        sasa_bb_desnat, sasa_sc_desnat = 100 .* gmx_sasa_per_restype(desnat_pdb, rname; gmx, n_dots, backbone, sidechain, ignore_hydrogen)
        sasas[rname] = OrderedDict(:sc => sasa_sc_desnat - sasa_sc_native, :bb => sasa_bb_desnat - sasa_bb_native)
    end
    return sasas # returns in Å^2
end

#
# Runs gmx sasa for a single residue type in a given PDB file.
#
function gmx_sasa_per_restype(
    pdbname, rname;
    gmx="gmx",
    n_dots::Int=500,
    backbone::Function,
    sidechain::Function,
    ignore_hydrogen::Bool=true,
)
    keepH(at) = ignore_hydrogen ? !(element(at) == "H") : true
    # Check if the gmx executable exists
    if isnothing(Sys.which(gmx))
        throw(ArgumentError("""\n
            Could not find GROMACS `$gmx` executable. Add the `gmx` command to the path or provide it explicitly with the `gmx` keyword argument.

        """))
    end
    index_file = tempname() * ".ndx"
    sasa_file = tempname() * ".xvg"
    p = read_pdb(pdbname, at -> isprotein(at) && keepH(at))
    inds_protein = index.(p)
    inds_sidechain = index.(select(p, at -> resname(at) == rname && sidechain(at)))
    inds_backbone = index.(select(p, at -> resname(at) == rname && backbone(at)))
    open(index_file, "w") do io
        write(io, "[ SC ]\n")
        for i in inds_sidechain
            write(io, "$i\n")
        end
        write(io, "[ BB ]\n")
        for i in inds_backbone
            write(io, "$i\n")
        end
        write(io, "[ PROT ]\n")
        for i in inds_protein
            write(io, "$i\n")
        end
    end
    if length(inds_sidechain) == 0 # special case for GLY
        sasa_sc = 0.0
        try
            run(pipeline(`$gmx sasa -s $pdbname -probe 0.14 -ndots $n_dots -surface PROT -output BB -n $index_file -o $sasa_file`, stdout=devnull, stderr=devnull))
        catch
            error("Error running $gmx sasa for of $resname in $pdbname")
        end
        sasa_bb = read_gmx_delta_sasa_per_restype_values(sasa_file, 1)[1]
    else
        try
            run(pipeline(`$gmx sasa -s $pdbname -probe 0.14 -ndots $n_dots -surface PROT -output SC BB -n $index_file -o $sasa_file`, stdout=devnull, stderr=devnull))
        catch
            error("Error running $gmx sasa for of $resname in $pdbname")
        end
        sasa_temp = read_gmx_delta_sasa_per_restype_values(sasa_file, 2)
        sasa_sc, sasa_bb = sasa_temp[1], sasa_temp[2]
    end
    rm(sasa_file, force=true)
    rm(index_file, force=true)
    return sasa_bb, sasa_sc
end

