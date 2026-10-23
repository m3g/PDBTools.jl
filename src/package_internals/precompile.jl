# Uses the test files distributed with PDBToolsBase (the other test files are lazy artifacts)
PrecompileTools.@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    PrecompileTools.@compile_workload begin
        ats = read_mmcif(PDBToolsBase.TESTCIF)
        ats = read_pdb(PDBToolsBase.SMALLPDB)
        select(ats, "protein")
        positions(ats)
        getseq(ats)
        center_of_mass(ats)
        eachresidue(ats)
        eachchain(ats)
        eachsegment(ats)
        eachmodel(ats)
        write_pdb(tempname(), ats)
        write_mmcif(tempname(), ats)
        residue_ticks(ats)
        moveto!(ats; center=[0.0, 0.0, 0.0])
        sasa(sasa_particles(ats))
    end
end