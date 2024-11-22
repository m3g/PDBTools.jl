PrecompileTools.@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    PrecompileTools.@compile_workload begin
        ats = read_mmcif(PDBTools.TESTCIF)
        ats = read_pdb(PDBTools.SMALLPDB)
        select(ats, "protein")
        coor(ats)
        getseq(ats)
        center_of_mass(ats)
        write_pdb(tempname(), ats)
        write_mmcif(tempname(), ats)
        residue_ticks(ats)
        moveto!(ats; center=[0.0, 0.0, 0.0])
    end
end