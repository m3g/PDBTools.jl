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
- `cosolvent::AbstractString`.
- `atoms::AbstractVector{<:PDBTools.Atom}`: Vector containing the atoms of the structure.
- `sasas::AbstractDict{String, AbstractDict{Symbol, Float64}}`: A dictionary containing the change in solvent accessible surface area (SASA)
  upon denaturation for each amino acid type. This data can be obtained from the `creamer_delta_sasa` function, the m-value server, or `delta_sasa_per_restype`:
    - The `creamer_delta_sasa` function provides estimated variations in SASA of a protein, using the Creamer unfolded model.
    - The output of the server can be parsed using the `parse_mvalue_server_sasa` function defined in this module.
    - Compute the SASA with `delta_sasa_per_restype`, a SASA calculation utility implemented in PDBTools.jl.
    - (deprecated) SASA values can be calculated using GROMACS with the `gmx_delta_sasa_per_restype` function.

- `type::Int`: Specifies which SASA value to use from the provided data, because the server provides minimum, average, and maximum values,
    according to different denatured models for the protein. The recommended value is `2` for comparison with experimental data.
    SASA values computed with `delta_sasa_per_restype` are single values, so `type=1` should be used in that case.

$(_available_cosolvents())

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
                parse_mvalue_server_sasa

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
```

## References

1. https://doi.org/10.1021/acs.jpcb.7b02138
2. https://doi.org/10.1016/s0076-6879(07)28023-1
3. https://www.pnas.org/doi/10.1073/pnas.0706251104

"""
function mvalue_delta_sasa(;
    model::Type{<:MValueModel}=MoeserHorinek,
    cosolvent::AbstractString="urea",
    atoms::AbstractVector{<:Atom}, sasas, type=1
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
- `unitcell=false`: `false`: no periodic boundary conditions; `true`: use the unit cell of the [`Structure`](@ref);
  or a [`UnitCell`](@ref).

"""
function delta_sasa_per_restype(;
    native::AbstractVector{<:Atom},
    desnat::AbstractVector{<:Atom},
    n_dots::Int=500,
    backbone::Function=at -> name(at) in ("N", "CA", "C", "O"),
    sidechain::Function=at -> !(name(at) in ("N", "CA", "C", "O")),
    ignore_hydrogen::Bool=true,
    unitcell::Union{Bool,UnitCell}=false,
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
