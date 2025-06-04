""" 
    select_with_vmd(atoms::AbstractVector{<:Atom}, selection::String; vmd="vmd", srcload=nothing)
    select_with_vmd(inputfile::String, selection::String; vmd="vmd", srcload=nothing)

Select atoms using vmd selection syntax, with vmd in background. The input can be a file or a list of atoms.

# Input structure and output format:

- `atoms::AbstractVector{<:Atom}`: A vector of `PDBTools.Atom` objects to select from. In this case, the
  output will be a vector of `PDBTools.Atom` objects that match the selection.
  
- `inputfile::String`: Path to the input file (e.g., PDB, PSF, GRO, etc.) or a temporary file containing atom data.
  In this case, two vectors will be returned: one with the indices of the selected atoms and another with their names.

*The outputs are different in each case because `VMD` supports selections on files like PSF, GRO, etc., which do not 
carry the full atom information like PDB files do.*

# Additional arguments:

- `selection::String`: A string containing the selection criteria in VMD syntax, e.g., `"protein and residue 1"`.
- `vmd::String`: The command to run VMD. Default is `"vmd"`, but can be set to the full path if VMD is not in the system PATH.
- `srcload::Union{Nothing, AbstractString, Vector{AbstractString}}`: A script or a list of VMD scripts to load before executing the selection,
  for example with macros to define custom selection keywords.

!!! warning
    VMD uses 0-based indexing, be careful to account for this when using indices in the selection string.
    If you want to suppress this warning, set `index_warning=false` in the function call.

"""
function select_with_vmd(inputfile::String, selection::String; vmd="vmd", srcload=nothing, index_warning=true)

    if index_warning
        if occursin("index", selection)
            @warn """\n

            Warning: VMD has 0-based indexing. 
            Selecting by index may not return atoms corresponding to the indices in the input file.

            To suppress this warning use `index_warning=false` in the function call.

            """
        end
    end

    if !isfile(inputfile)
        error("Could not find file: $inputfile")
    end

    vmdinput_file = tempname()
    vmd_input = Base.open(vmdinput_file, "w")
    if !isnothing(srcload)
        if srcload isa AbstractString
            srcload = [srcload]
        end
        for srcfile in srcload
            Base.write(vmd_input, "source \"$srcfile\" \n")
        end
    end
    Base.write(vmd_input, "mol new \"$inputfile\" \n")
    Base.write(vmd_input, "set sel [ atomselect top \"$selection\" ] \n")
    Base.write(vmd_input, "puts \"INDEXLIST\" \n")
    Base.write(vmd_input, "set indices [ \$sel get index ] \n")
    Base.write(vmd_input, "puts \"ENDINDEXLIST\" \n")
    Base.write(vmd_input, "puts \"NAMELIST\" \n")
    Base.write(vmd_input, "set names [ \$sel get name ] \n")
    Base.write(vmd_input, "puts \"ENDNAMELIST\" \n")
    Base.write(vmd_input, "exit \n")
    Base.close(vmd_input)

    vmd_output = Base.read(`$vmd -dispdev text -e $vmdinput_file`, String)

    # Read indices
    local index_list::String
    readnext = false
    for line in split(vmd_output, "\n")
        if occursin("cannot parse selection text", line)
            error("cannot parse selection: \"$selection\".")
        end
        if readnext
            if line == "ENDINDEXLIST"
                error("selection \"$selection\" is empty.")
            end
            index_list = line
            break
        end
        if line == "INDEXLIST"
            readnext = true
        end
    end
    index_split = split(index_list)
    nsel = length(index_split)
    selection_indices = Vector{Int}(undef, nsel)
    for i = 1:nsel
        selection_indices[i] = parse(Int, index_split[i]) + 1
    end

    # Read atom names
    local name_list::String
    readnext = false
    for line in split(vmd_output, "\n")
        if readnext
            if line == "ENDNAMELIST"
                error("ERROR: Selection '$selection' does not contain any atom")
            end
            name_list = line
            break
        end
        if line == "NAMELIST"
            readnext = true
        end
    end
    name_split = split(name_list)
    nsel = length(name_split)
    selection_names = Vector{String}(undef, nsel)
    for i = 1:nsel
        selection_names[i] = strip(name_split[i])
    end

    return selection_indices, selection_names
end

function select_with_vmd(atoms::AbstractVector{<:Atom}, selection::String; vmd="vmd", srcload=nothing)
    tmp_file = tempname()
    write_pdb(tmp_file, atoms)
    inds, _ = select_with_vmd(tmp_file, selection; vmd=vmd, srcload=srcload)
    return select(atoms, at -> index(at) in inds)
end

@testitem "select_with_vmd" begin
    pdbfile = PDBTools.TESTPDB
    if !isnothing(Sys.which("vmd"))
        inds, names = select_with_vmd(pdbfile, "protein and residue 1")
        @test inds == [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
        @test names == ["N", "HN", "CA", "HA", "CB", "HB1", "HB2", "SG", "HG1", "C", "O"]
        atoms = read_pdb(pdbfile)
        sel = select_with_vmd(atoms, "protein and residue 1")
        @test index.(sel) == [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
        @test name.(sel) == ["N", "HN", "CA", "HA", "CB", "HB1", "HB2", "SG", "HG1", "C", "O"]
    end
end
