@testitem "mvalue" begin
    using PDBTools
    using ShowMethodTesting

    dir = @__DIR__
    MJC_native = read_pdb(joinpath(dir, "1MJC_native.pdb"), "protein")
    MJC_desnat = read_pdb(joinpath(dir, "1MJC_straight.pdb"), "protein")

    r_1MJC = mvalue(MJC_native, MJC_desnat, "urea")
    @test isapprox(r_1MJC.tot, -1.334; atol=1e-2)
    @test isapprox(r_1MJC.bb, -1.405; atol=1e-2)
    @test isapprox(r_1MJC.sc, 0.071; atol=1e-2)

    # Confirm that now cosolvent selection is not case-sensitive
    r_1MJC = mvalue(MJC_native, MJC_desnat, "Urea")
    @test isapprox(r_1MJC.tot, -1.334; atol=1e-2)

    # test non-parallel vs. parallel calculations
    ms = mvalue(MJC_native, MJC_desnat, "urea"; parallel=false)
    mp = mvalue(MJC_native, MJC_desnat, "urea"; parallel=true)
    @test ms.tot ≈ mp.tot

    # Test show method
    @test parse_show(r_1MJC; repl=Dict(r"PDBTools." => "")) ≈ """
        MValue{AutonBolen} - 69 residues - cosolvent: "urea"
            Total m-value: -1.3340737 kcal mol⁻¹
            Backbone contributions: -1.4051082 kcal mol⁻¹
            Side-chain contributions: 0.07103452 kcal mol⁻¹
        """

    # Test save/load preserves model type and data
    tmp = tempname() * ".json"
    save(tmp, r_1MJC)
    m_load = load(MValue, tmp)
    @test m_load isa MValue{AutonBolen}
    @test m_load.nresidues == r_1MJC.nresidues
    @test m_load.tot == r_1MJC.tot
    @test m_load.bb == r_1MJC.bb
    @test m_load.sc == r_1MJC.sc
    @test m_load.residue_contributions_bb == r_1MJC.residue_contributions_bb
    @test m_load.residue_contributions_sc == r_1MJC.residue_contributions_sc
    @test m_load.cosolvent == r_1MJC.cosolvent
    rm(tmp; force=true)

    r_1MJC = mvalue(MJC_native, MJC_desnat, "tmao")
    @test isapprox(r_1MJC.tot, 2.5198; atol=1e-2)
    @test isapprox(r_1MJC.bb, 3.2425; atol=1e-2)
    @test isapprox(r_1MJC.sc, -0.7227; atol=1e-2)

    #
    # MoeserHorinek model
    # 
    r_1MJC = mvalue(MJC_native, MJC_desnat, "urea"; model=MoeserHorinek)
    @test isapprox(r_1MJC.tot, -1.203; atol=1e-2)
    @test isapprox(r_1MJC.bb, -0.676; atol=1e-2)
    @test isapprox(r_1MJC.sc, -0.526; atol=1e-2)

    save(tmp, r_1MJC)
    m_load = load(MValue, tmp)
    @test m_load isa MValue{MoeserHorinek}
    @test m_load.tot == r_1MJC.tot
    rm(tmp; force=true)

    # Same definitions as Gromacs for side chain and backbone
    r_1MJC = mvalue(MJC_native, MJC_desnat, "urea";
        model=MoeserHorinek,
        backbone=at -> name(at) in ("N", "CA", "C", "O"),
        sidechain=at -> !(name(at) in ("N", "CA", "C", "O")),
    )
    @test isapprox(r_1MJC.tot, -1.2025; atol=1e-2)
    @test isapprox(r_1MJC.bb, -0.676; atol=1e-2)
    @test isapprox(r_1MJC.sc, -0.526; atol=1e-2)

    # Compute using precomputed SASAs
    sn = sasa_particles(CreamerUnitedAtomRadii, MJC_native)
    sd = sasa_particles(CreamerUnitedAtomRadii, MJC_desnat)
    r_1MJC = mvalue(sn, sd, "urea"; model=MoeserHorinek)
    @test isapprox(r_1MJC.tot, -1.2025; atol=1e-2)
    @test isapprox(r_1MJC.bb, -0.6759; atol=1e-2)
    @test isapprox(r_1MJC.sc, -0.52656; atol=1e-2)

    # Provide a selection
    r_1MJC = mvalue(MJC_native, MJC_desnat, "urea"; sel="acidic")
    @test isapprox(r_1MJC.tot, -0.01512; atol=1e-2)
    @test isapprox(r_1MJC.bb, -0.10998993; atol=1e-2)
    @test isapprox(r_1MJC.sc, 0.09485909; atol=1e-2)

    # Test selection additivity
    r0 = mvalue(MJC_native, MJC_desnat, "urea")
    r1 = mvalue(MJC_native, MJC_desnat, "urea"; sel="acidic")
    r2 = mvalue(MJC_native, MJC_desnat, "urea"; sel="not acidic")
    @test r0.nresidues == r1.nresidues + r2.nresidues
    @test r0.tot ≈ r1.tot + r2.tot
    rbb = mvalue(MJC_native, MJC_desnat, "urea"; sel="backbone")
    @test r0.bb ≈ rbb.tot
    @test r0.bb ≈ rbb.bb
    @test rbb.sc ≈ 0.0 atol=1e-10
    rsc = mvalue(MJC_native, MJC_desnat, "urea"; sel="sidechain")
    @test r0.sc ≈ rsc.tot
    @test r0.sc ≈ rsc.sc
    @test rsc.bb ≈ 0.0 atol=1e-10

    # Provide a selection using SASAs
    r_1MJC = mvalue(sn, sd, "urea"; sel="acidic")
    @test isapprox(r_1MJC.tot, -0.01512; atol=1e-2)
    @test isapprox(r_1MJC.bb, -0.10998993; atol=1e-2)
    @test isapprox(r_1MJC.sc, 0.09486909; atol=1e-2)

    # Test selection additivity
    r1 = mvalue(sn, sd, "urea"; sel="acidic")
    r2 = mvalue(sn, sd, "urea"; sel="not acidic")
    @test r0.nresidues == r1.nresidues + r2.nresidues
    @test r0.tot ≈ r1.tot + r2.tot

    # Some errors
    p1 = read_pdb(PDBTools.TESTPDB, "protein and residue < 10")
    p2 = read_pdb(PDBTools.TESTPDB, "protein and residue < 11")
    @test_throws "same number of residues" mvalue(p1, p2, "urea") # different number of residues

    p1 = select(p1, "protein")
    p2 = copy.(p1)
    rs = collect(eachresidue(p2))
    for at in rs[7]
        at.resname = "ASP"
    end
    @test_throws "same type" mvalue(p1, p2, "urea") # different number of residues

    # Support for periodic boundary conditions
    no_pbc = read_pdb(PDBTools.TESTNOPBC, "protein")
    pbc = read_pdb(PDBTools.TESTPBC, "protein")
    uc = [107.845, 107.845, 107.845]
    @test transfer_free_energy(no_pbc, "urea").tot ≈ 
          transfer_free_energy(pbc, "urea"; unitcell=uc).tot
    @test mvalue(pbc, pbc, "urea"; unitcell=uc).tot ≈ 0.0 atol=1e-5

    # Error if wrong atomic radii was provided
    s = sasa_particles(no_pbc)
    @test_throws "must use CreamerUnitedAtomRadii" mvalue(s, s, "urea")

end

@testitem "mvalue_delta_sasa" begin
    using PDBTools: mvalue_delta_sasa
    using PDBTools: delta_sasa_per_restype
    using PDBTools: parse_mvalue_server_sasa
    using PDBTools: gmx_delta_sasa_per_restype

    sasa_1MJC_clean = parse_mvalue_server_sasa(
        """
        ALA 		    5 	 (    95.6)    138.1 [   180.6] 	 | 	 (    -4.6)     35.7 [    75.9] 
        PHE 		    6 	 (   430.8)    491.1 [   557.4] 	 | 	 (    65.0)    119.0 [   173.0] 
        LEU 		    2 	 (   131.2)    147.5 [   163.8] 	 | 	 (   -29.8)    -13.8 [     2.2] 
        ILE 		    4 	 (   247.0)    315.0 [   383.0] 	 | 	 (    24.9)     43.9 [    62.9] 
        VAL 		    5 	 (   403.1)    475.9 [   548.6] 	 | 	 (    76.9)     99.4 [   121.9] 
        PRO 		    2 	 (   134.0)    141.0 [   148.0] 	 | 	 (    27.2)     34.4 [    41.6] 
        MET 		    1 	 (    55.6)     72.7 [    89.8] 	 | 	 (    16.1)     24.6 [    33.2] 
        TRP 		    1 	 (    71.5)     73.3 [    75.2] 	 | 	 (    14.5)     22.9 [    31.4] 
        GLY 		   10 	 (     0.0)      0.0 [     0.0] 	 | 	 (   200.9)    306.4 [   411.9] 
        SER 		    7 	 (   110.3)    157.9 [   205.5] 	 | 	 (   -49.0)     -9.8 [    29.4] 
        THR 		    4 	 (   130.9)    158.7 [   186.5] 	 | 	 (    70.0)     91.8 [   113.6] 
        TYR 		    1 	 (    76.4)     87.0 [    97.7] 	 | 	 (    -2.1)      5.8 [    13.7] 
        GLN 		    2 	 (    89.8)    113.5 [   137.2] 	 | 	 (    18.8)     35.0 [    51.2] 
        ASN 		    3 	 (    54.5)     71.1 [    87.8] 	 | 	 (    43.2)     65.9 [    88.5] 
        ASP 		    6 	 (    70.4)    117.2 [   164.0] 	 | 	 (   -53.0)     -5.6 [    41.8] 
        GLU 		    2 	 (    30.2)     51.3 [    72.4] 	 | 	 (    -4.5)     11.1 [    26.7] 
        HIS 		    1 	 (    42.7)     50.3 [    57.9] 	 | 	 (    13.2)     22.4 [    31.7] 
        LYS 		    7 	 (   260.6)    317.6 [   374.7] 	 | 	 (   -10.0)     44.3 [    98.5] 
        ARG 		    0 	 (     0.0)      0.0 [     0.0] 	 | 	 (     0.0)      0.0 [     0.0] 
        CYS 		    0 	 (     0.0)      0.0 [     0.0] 	 | 	 (     0.0)      0.0 [     0.0] 
        """
    )

    sasa_2RN2_clean = parse_mvalue_server_sasa(
        """
        ALA 		   14 	 (   426.2)    545.2 [   664.2] 	 | 	 (   149.9)    262.6 [   375.3] 
        PHE 		    2 	 (   187.5)    207.6 [   229.7] 	 | 	 (    29.9)     47.9 [    65.9] 
        LEU 		   12 	 (   929.1)   1026.9 [  1124.7] 	 | 	 (    91.9)    187.9 [   283.9] 
        ILE 		    7 	 (   603.9)    722.9 [   841.9] 	 | 	 (    55.7)     89.0 [   122.2] 
        VAL 		    9 	 (   518.5)    649.4 [   780.4] 	 | 	 (    76.4)    116.9 [   157.4] 
        PRO 		    5 	 (   161.5)    179.0 [   196.5] 	 | 	 (    51.9)     69.9 [    87.9] 
        MET 		    4 	 (   187.8)    256.2 [   324.6] 	 | 	 (   -29.5)      4.7 [    38.9] 
        TRP 		    6 	 (   716.2)    727.3 [   738.4] 	 | 	 (    29.7)     80.4 [   131.1] 
        GLY 		   14 	 (     0.0)      0.0 [     0.0] 	 | 	 (   439.1)    586.8 [   734.5] 
        SER 		    4 	 (   197.8)    225.0 [   252.2] 	 | 	 (    54.6)     77.0 [    99.4] 
        THR 		   10 	 (   397.0)    466.5 [   536.0] 	 | 	 (    37.0)     91.5 [   146.0] 
        TYR 		    5 	 (   556.9)    610.1 [   663.4] 	 | 	 (    65.4)    104.9 [   144.4] 
        GLN 		    8 	 (   143.6)    238.4 [   333.2] 	 | 	 (    67.8)    132.6 [   197.4] 
        ASN 		    7 	 (   182.6)    221.5 [   260.3] 	 | 	 (    76.7)    129.6 [   182.4] 
        ASP 		    7 	 (   306.3)    360.9 [   415.5] 	 | 	 (    48.7)    104.0 [   159.3] 
        GLU 		   12 	 (   426.2)    552.8 [   679.4] 	 | 	 (    90.4)    184.0 [   277.6] 
        HIS 		    5 	 (    87.1)    125.1 [   163.1] 	 | 	 (   -20.2)     26.0 [    72.3] 
        LYS 		   11 	 (   317.2)    406.8 [   496.5] 	 | 	 (    48.7)    134.0 [   219.2] 
        ARG 		   10 	 (   546.2)    688.2 [   830.2] 	 | 	 (   110.1)    189.6 [   269.1] 
        CYS 		    3 	 (   183.2)    213.4 [   243.5] 	 | 	 (    32.3)     56.7 [    81.2] 
        """
    )

    # Result from Moeser and Horinek https://doi.org/10.1021/jp409934q
    references = Dict(
        #                    total                  bb                   sc         
        "1MJC" => (-0.8470588235294114, -0.44117647058823506, -0.44117647058823506),
        "2RN2" => (-2.3647058823529408, -1.1999999999999997, -1.164705882352941)
    )

    gmx = Sys.which("gmx")
    if isnothing(gmx)
        @warn "gmx executable not available: some tests won't be run"
    end

    dir = @__DIR__

    #
    # 1MJC
    #
    MJC_native = read_pdb(joinpath(dir, "1MJC_native.pdb"), "protein")
    MJC_desnat = read_pdb(joinpath(dir, "1MJC_straight.pdb"), "protein")

    r_1MJC = mvalue_delta_sasa(; atoms=MJC_native, sasas=sasa_1MJC_clean, type=2)
    @test isapprox(r_1MJC.tot, references["1MJC"][1]; rtol=1e-1)
    @test isapprox(r_1MJC.bb, references["1MJC"][2]; rtol=1e-1)
    @test isapprox(r_1MJC.sc, references["1MJC"][3]; rtol=1e-1)

    # with PDBTools.sasa
    sasa_1MJC_julia = delta_sasa_per_restype(; native=MJC_native, desnat=MJC_desnat)
    r_1MJC = mvalue_delta_sasa(; atoms=MJC_native, sasas=sasa_1MJC_julia)
    @test isapprox(r_1MJC.tot, -1.22; rtol=1e-1)
    @test isapprox(r_1MJC.bb, -0.755; rtol=1e-1)
    @test isapprox(r_1MJC.sc, -0.468; rtol=1e-1)
    sasa_1MJC_julia = delta_sasa_per_restype(; native=MJC_native, desnat=MJC_desnat, ignore_hydrogen=false)
    r_1MJC = mvalue_delta_sasa(; atoms=MJC_native, sasas=sasa_1MJC_julia)
    @test isapprox(r_1MJC.tot, -0.9863; rtol=1e-1)
    @test isapprox(r_1MJC.bb, -0.344; rtol=1e-1)
    @test isapprox(r_1MJC.sc, -0.642; rtol=1e-1)

    if !isnothing(gmx)
        sasa_1MJC = gmx_delta_sasa_per_restype(;
            native_pdb=joinpath(dir, "1MJC_native.pdb"),
            desnat_pdb=joinpath(dir, "1MJC_straight.pdb"),
        )
        r_1MJC = mvalue_delta_sasa(; atoms=MJC_native, sasas=sasa_1MJC)
        @test isapprox(r_1MJC.tot, -1.2335; rtol=1e-1)
        @test isapprox(r_1MJC.bb, -0.7333; rtol=1e-1)
        @test isapprox(r_1MJC.sc, -0.50025; rtol=1e-1)
        sasa_1MJC = gmx_delta_sasa_per_restype(;
            native_pdb=joinpath(dir, "1MJC_native.pdb"),
            desnat_pdb=joinpath(dir, "1MJC_straight.pdb"),
            ignore_hydrogen=false,
        )
        r_1MJC = mvalue_delta_sasa(; atoms=MJC_native, sasas=sasa_1MJC)
        @test isapprox(r_1MJC.tot, -0.98749; rtol=1e-1)
        @test isapprox(r_1MJC.bb, -0.3014; rtol=1e-1)
        @test isapprox(r_1MJC.sc, -0.67707; rtol=1e-1)
    end

    #
    # 2RN2
    #
    RN2_native = read_pdb(joinpath(dir, "2RN2_native.pdb"), "protein")
    RN2_desnat = read_pdb(joinpath(dir, "2RN2_straight.pdb"), "protein")

    r_2RN2 = mvalue_delta_sasa(; atoms=RN2_native, sasas=sasa_2RN2_clean, type=2)
    @test isapprox(r_2RN2.tot, references["2RN2"][1]; rtol=1e-1)
    @test isapprox(r_2RN2.bb, references["2RN2"][2]; rtol=1e-1)
    @test isapprox(r_2RN2.sc, references["2RN2"][3]; rtol=1e-1)

    sasa_2RN2_julia = delta_sasa_per_restype(; native=RN2_native, desnat=RN2_desnat)
    r_2RN2 = mvalue_delta_sasa(; atoms=RN2_native, sasas=sasa_2RN2_julia)
    @test isapprox(r_2RN2.tot, -3.13; rtol=1e-1)
    @test isapprox(r_2RN2.bb, -1.89; rtol=1e-1)
    @test isapprox(r_2RN2.sc, -1.23; rtol=1e-1)

    sasa_2RN2_julia = delta_sasa_per_restype(; native=RN2_native, desnat=RN2_desnat, ignore_hydrogen=false)
    r_2RN2 = mvalue_delta_sasa(; atoms=RN2_native, sasas=sasa_2RN2_julia)
    @test isapprox(r_2RN2.tot, -2.59; rtol=1e-1)
    @test isapprox(r_2RN2.bb, -1.03; rtol=1e-1)
    @test isapprox(r_2RN2.sc, -1.57; rtol=1e-1)

    if !isnothing(gmx)
        sasa_2RN2 = gmx_delta_sasa_per_restype(;
            native_pdb=joinpath(dir, "2RN2_native.pdb"),
            desnat_pdb=joinpath(dir, "2RN2_straight.pdb"),
        )
        r_2RN2 = mvalue_delta_sasa(; atoms=RN2_native, sasas=sasa_2RN2)
        @test isapprox(r_2RN2.tot, -3.17; rtol=1e-1)
        @test isapprox(r_2RN2.bb, -1.94; rtol=1e-1)
        @test isapprox(r_2RN2.sc, -1.23; rtol=1e-1)
        sasa_2RN2 = gmx_delta_sasa_per_restype(;
            native_pdb=joinpath(dir, "2RN2_native.pdb"),
            desnat_pdb=joinpath(dir, "2RN2_straight.pdb"),
            ignore_hydrogen=false,
        )
        r_2RN2 = mvalue_delta_sasa(; atoms=RN2_native, sasas=sasa_2RN2)
        @test isapprox(r_2RN2.tot, -2.54; rtol=1e-1)
        @test isapprox(r_2RN2.bb, -0.88; rtol=1e-1)
        @test isapprox(r_2RN2.sc, -1.66; rtol=1e-1)
    end

    # Results from running the Auton and Bolen m-value server with 1MJC_clean.pdb
    # http://best.bio.jhu.edu/mvalue/
    data_mvalue_server_auton_and_bolen_1MJC = Dict(
        "tmao" => (224, 1160, 2095),
        "sarcosine" => (406, 1010, 1614),
        "betaine" => (-502, 76, 650),
        "proline" => (-200, 226, 649),
        "sorbitol" => (378, 780, 1183),
        "sucrose" => (189, 876, 1559),
        "urea" => (-293, -711, -1132),
    )

    MJC_clean = read_pdb(joinpath(dir, "1MJC_clean.pdb"))
    for (cos, dg) in data_mvalue_server_auton_and_bolen_1MJC
        for ig in 1:3
            m = mvalue_delta_sasa(model=AutonBolen, cosolvent=cos, atoms=MJC_clean, sasas=sasa_1MJC_clean, type=ig)
            @test m.tot ≈ 1e-3 * dg[ig] rtol = 0.1
            cm = CreamerDenaturedModel(MJC_clean, ig)
            m = mvalue(cm, cos)
            @test m.tot ≈ 1e-3 * dg[ig] atol = 0.1
            @test sum(m.residue_contributions_bb) + sum(m.residue_contributions_sc) ≈ 1e-3 * dg[ig] atol = 0.1
        end
    end

    data_mvalue_server_auton_and_bolen_2RN2 = Dict(
        "tmao" => (978, 3215, 5453),
        "sarcosine" => (1410, 2867, 4323),
        "betaine" => (-1301, 130, 1560),
        "proline" => (-568, 464, 1496),
        "sorbitol" => (830, 1757, 2684),
        "sucrose" => (952, 2578, 4202),
        "urea" => (-1217, -2226, -3236)
    )

    RN2_clean = read_pdb(joinpath(dir, "2RN2_clean.pdb"))
    for (cos, dg) in data_mvalue_server_auton_and_bolen_2RN2
        for ig in 1:3
            m = mvalue_delta_sasa(model=AutonBolen, cosolvent=cos, atoms=RN2_clean, sasas=sasa_2RN2_clean, type=ig)
            @test m.tot ≈ 1e-3 * dg[ig] rtol = 0.1
            cm = CreamerDenaturedModel(RN2_clean, ig)
            m = mvalue(cm, cos)
            @test m.tot ≈ 1e-3 * dg[ig] atol = 0.1
        end
    end

    # gmx keyword argument tests
    @test_throws ArgumentError gmx_delta_sasa_per_restype(;
        native_pdb=joinpath(dir, "1MJC_native.pdb"),
        desnat_pdb=joinpath(dir, "1MJC_straight.pdb"),
        gmx="/tmp/gmx_fake",
    )
    if !isnothing(gmx)
        @test length(gmx_delta_sasa_per_restype(;
            native_pdb=joinpath(dir, "2RN2_native.pdb"),
            desnat_pdb=joinpath(dir, "2RN2_straight.pdb"),
            gmx=gmx,
        )) == 20
    end

    # Test m-value calculation using PBCs
    pbc = read_pdb(joinpath(dir, "pbc.pdb"))
    no_pbc = read_pdb(joinpath(dir, "no_pbc.pdb"))
    uc = lattice_to_matrix(107.845413, 107.845413, 107.845413, 90., 90., 90.)
    ss = delta_sasa_per_restype(; native=pbc, desnat=no_pbc, unitcell=uc)
    m = mvalue_delta_sasa(model=AutonBolen, cosolvent="urea", atoms=pbc, sasas=ss)
    @test m.tot ≈ 0.0 atol = 1e-3

    #=

    Average SASA values from lower and upper bounds of the denatured state ensemble,
    reported in Supplementary Table 2 of https://doi.org/10.1073/pnas.0507053102
    These are the average values of Table 1 of https://doi.org/10.1021/bi962819o
    Can be used for testing, but **do not** provide the same accuracy as the values
    calculated with GROMACS or obtained from the m-value server.

    =#
    const sasa_desnat_average = Dict(
        "ALA" => Dict(:bb => 27.9, :sc => 55.1),
        "PHE" => Dict(:bb => 24.3, :sc => 128.8),
        "LEU" => Dict(:bb => 22.7, :sc => 109.6),
        "ILE" => Dict(:bb => 20.0, :sc => 117.1),
        "VAL" => Dict(:bb => 20.4, :sc => 96.4),
        "PRO" => Dict(:bb => 22.5, :sc => 87.0),
        "MET" => Dict(:bb => 25.3, :sc => 122.4),
        "TRP" => Dict(:bb => 23.6, :sc => 156.6),
        "GLY" => Dict(:bb => 65.2, :sc => 0.0),
        "SER" => Dict(:bb => 29.4, :sc => 66.5),
        "THR" => Dict(:bb => 24.1, :sc => 84.3),
        "TYR" => Dict(:bb => 25.6, :sc => 141.7),
        "GLN" => Dict(:bb => 25.3, :sc => 116.9),
        "ASN" => Dict(:bb => 25.2, :sc => 90.1),
        "ASP" => Dict(:bb => 26.0, :sc => 87.0),
        "GLU" => Dict(:bb => 25.7, :sc => 113.4),
        "HIS" => Dict(:bb => 24.2, :sc => 111.5),
        "LYS" => Dict(:bb => 26.1, :sc => 150.7),
        "ARG" => Dict(:bb => 25.1, :sc => 171.1),
        "CYS" => Dict(:bb => 26.4, :sc => 73.0),
    )

end

@testitem "creamer_delta_sasa" begin
    using PDBTools
    using PDBTools: mvalue_delta_sasa, creamer_delta_sasa

    dir = @__DIR__

    # Results from running the Auton and Bolen m-value server with 1MJC_clean.pdb
    data_mvalue_server_auton_and_bolen_1MJC = Dict(
        "tmao" => (224, 1160, 2095),
        "sarcosine" => (406, 1010, 1614),
        "betaine" => (-502, 76, 650),
        "proline" => (-200, 226, 649),
        "sorbitol" => (378, 780, 1183),
        "sucrose" => (189, 876, 1559),
        "urea" => (-293, -711, -1132),
    )

    MJC_clean = read_pdb(joinpath(dir, "1MJC_clean.pdb"))
    for (cos, dg) in data_mvalue_server_auton_and_bolen_1MJC
        for ig in 1:3
            m = mvalue_delta_sasa(
                model=AutonBolen,
                cosolvent=cos,
                atoms=MJC_clean,
                sasas=creamer_delta_sasa(MJC_clean),
                type=ig
            )
            @test m.tot ≈ 1e-3 * dg[ig] atol = 0.1
        end
    end

    data_mvalue_server_auton_and_bolen_2RN2 = Dict(
        "tmao" => (978, 3215, 5453),
        "sarcosine" => (1410, 2867, 4323),
        "betaine" => (-1301, 130, 1560),
        "proline" => (-568, 464, 1496),
        "sorbitol" => (830, 1757, 2684),
        "sucrose" => (952, 2578, 4202),
        "urea" => (-1217, -2226, -3236)
    )

    RN2_clean = read_pdb(joinpath(dir, "2RN2_clean.pdb"))
    for (cos, dg) in data_mvalue_server_auton_and_bolen_2RN2
        for ig in 1:3
            m = mvalue_delta_sasa(
                model=AutonBolen,
                cosolvent=cos,
                atoms=RN2_clean,
                sasas=creamer_delta_sasa(RN2_clean),
                type=ig
            )
            @test m.tot ≈ 1e-3 * dg[ig] atol = 0.1
        end
    end

    # test that hydrogens are properly handled (0 vdw radius)
    pdb = read_pdb(PDBTools.TESTPDB, "protein or name CLA")
    prot = select(pdb, "protein")
    mH = mvalue_delta_sasa(; atoms=prot, sasas=PDBTools.creamer_delta_sasa(prot))
    prot = select(pdb, "protein and not element H")
    m_not_H = mvalue_delta_sasa(; atoms=prot, sasas=PDBTools.creamer_delta_sasa(prot))
    @test mH.tot ≈ m_not_H.tot

    # test error for non-recognized element
    @test_throws "Could not determine" PDBTools.creamer_delta_sasa(pdb)

end

@testitem "transfer free energy" begin
    using PDBTools
    using ShowMethodTesting
    dir = @__DIR__

    #
    # Server results for 2RN2_clean.pdb
    #
    #     Native State Transfer Free Energy Contributions to the m-value for /tmp/phpt5bGI3. 
    server_result = Dict(
        "tmao" => (bb=3573, sc=-1945, tot=1629),
        "sarcosine" => (bb=2065, sc=-715, tot=1349),
        "betaine" => (bb=2660, sc=-3334, tot=-674),
        "proline" => (bb=1906, sc=-2207, tot=-301),
        "Glycerol" => (bb=556, sc=-1093, tot=-537),
        "sorbitol" => (bb=1390, sc=-1011, tot=379),
        "sucrose" => (bb=2462, sc=-1698, tot=764),
        "Trehalose" => (bb=2462, sc=-1341, tot=1121),
        "urea" => (bb=-1548, sc=445, tot=-1104),
    )
    p = read_pdb(joinpath(dir, "2RN2_clean.pdb"))
    for (cosolvent, vals)  in pairs(server_result)
        t = transfer_free_energy(p, cosolvent)
        @test t.bb ≈ 1e-3*vals.bb atol=0.1
        @test t.sc ≈ 1e-3*vals.sc atol=0.1
        @test t.tot ≈ 1e-3*vals.tot atol=0.1
        # from SASAs
        s = sasa_particles(CreamerUnitedAtomRadii, p)
        t = transfer_free_energy(s, cosolvent)
        @test t.bb ≈ 1e-3*vals.bb atol=0.1
        @test t.sc ≈ 1e-3*vals.sc atol=0.1
        @test t.tot ≈ 1e-3*vals.tot atol=0.1
    end

    # Test selection additivity
    r0 = transfer_free_energy(p, "urea")
    r1 = transfer_free_energy(p, "urea"; sel="polar")
    r2 = transfer_free_energy(p, "urea"; sel="not polar")
    @test r0.nresidues == r1.nresidues + r2.nresidues
    @test r0.tot ≈ r1.tot + r2.tot

    # Provide a selection using SASAs
    s = sasa_particles(CreamerUnitedAtomRadii, p)
    r0 = transfer_free_energy(s, "urea")
    r1 = transfer_free_energy(s, "urea"; sel="polar")
    r2 = transfer_free_energy(s, "urea"; sel="not polar")
    @test r0.nresidues == r1.nresidues + r2.nresidues
    @test r0.tot ≈ r1.tot + r2.tot

    # Test selection additivity
    r1 = transfer_free_energy(s, "urea"; sel="polar")
    r2 = transfer_free_energy(s, "urea"; sel="not polar")
    @test r0.nresidues == r1.nresidues + r2.nresidues
    @test r0.tot ≈ r1.tot + r2.tot
    rbb = transfer_free_energy(s, "urea"; sel="backbone")
    @test r0.bb ≈ rbb.tot
    @test r0.bb ≈ rbb.bb
    @test rbb.sc ≈ 0.0 atol=1e-10
    rsc = transfer_free_energy(s, "urea"; sel="sidechain")
    @test r0.sc ≈ rsc.tot
    @test r0.sc ≈ rsc.sc
    @test rsc.bb ≈ 0.0 atol=1e-10

    # Test show method
    t = transfer_free_energy(p, "urea")
    @test parse_show(t; repl=Dict(r"PDBTools." => "")) ≈ """
        PDBTools.TransferFreeEnergy{AutonBolen} - 155 residues to 1M "urea".
        Total transfer free energy: -1.125963 kcal mol⁻¹
        Backbone contributions: -1.5711304 kcal mol⁻¹
        Side-chain contributions: 0.4451674 kcal mol⁻¹
        """

    # Test save/load preserves model type and data
    tmp = tempname() * ".json"
    save(tmp, t)
    t_load = load(TransferFreeEnergy, tmp)
    @test t_load isa TransferFreeEnergy{AutonBolen}
    @test t_load.nresidues == t.nresidues
    @test t_load.tot == t.tot
    @test t_load.bb == t.bb
    @test t_load.sc == t.sc
    @test t_load.residue_contributions_bb == t.residue_contributions_bb
    @test t_load.residue_contributions_sc == t.residue_contributions_sc
    @test t_load.cosolvent == t.cosolvent
    rm(tmp; force=true)

    t_mh = transfer_free_energy(p, "sucrose"; model=MoeserHorinek)
    save(tmp, t_mh)
    t_mh_load = load(TransferFreeEnergy, tmp)
    @test t_mh_load isa TransferFreeEnergy{MoeserHorinek}
    @test t_mh_load.tot == t_mh.tot
    rm(tmp; force=true)

    # Tests for MoeserHorinekApp model
    cm = CreamerDenaturedModel(read_pdb(PDBTools.TESTPDB, "protein"))
    for c in filter(!=("urea"), keys(PDBTools.cosolvent_column_MoeserHorinekApp))
        c_ab = mvalue(cm, c; model=AutonBolen)
        c_mhapp = mvalue(cm, c; model=MoeserHorinekApp)
        @test c_ab.tot ≈ c_mhapp.tot atol=0.25
    end

    # Tests for Accessibility model
    cm = CreamerDenaturedModel(read_pdb(PDBTools.TESTPDB, "protein"))
    for c in keys(PDBTools.cosolvent_column_Accessibility)
        c_ab = mvalue(cm, c; model=AutonBolen).tot
        c_acc = mvalue(cm, c; model=Accessibility).tot
        @test c_ab ≈ c_acc atol = 0.25
        @test 1 + c_ab ≈ 1 + c_acc rtol = 0.25
    end

    # Tests for MTRecord model (native structure vs. its own fully-extended chain).
    # For urea, Guinn et al. (Table S3) themselves use this exact "extended beta (all-trans)"
    # reference state, so these reproduce the paper's own predicted m-values well.
    RN2 = read_pdb(joinpath(dir, "2RN2_native.pdb"), "protein")
    MJC = read_pdb(joinpath(dir, "1MJC_clean.pdb"))
    rec_m = MTRecordDenaturedModel(RN2)
    c_rec = mvalue(rec_m, "urea").tot
    @test c_rec ≈ -2.2 rtol = 0.15
    rec_m = MTRecordDenaturedModel(MJC)
    c_rec = mvalue(rec_m, "urea").tot
    @test c_rec ≈ -0.94 rtol = 0.1

    # Test providing models or chains to MTRecordDenaturedModel
    MJC_m = collect(eachmodel(MJC))[1]
    MJC_c = collect(eachchain(MJC))[1]
    @test mvalue(MTRecordDenaturedModel(MJC_m), "urea").tot ≈ c_rec
    @test mvalue(MTRecordDenaturedModel(MJC_c), "urea").tot ≈ c_rec

    # Test consistency of the TFE path: mvalue(MTRecordDenaturedModel(MJC), ...) must agree
    # with manually differencing the transfer free energies of MJC and its extended chain.
    c_rec_mjc = mvalue(MTRecordDenaturedModel(MJC), "urea").tot
    tfe_n = transfer_free_energy(MJC, "urea"; model=MTRecord)
    tfe_d = transfer_free_energy(extended_chain(MJC), "urea"; model=MTRecord)
    m_tfe = tfe_d.tot - tfe_n.tot
    @test c_rec_mjc ≈ m_tfe

    # Test path of TFE from SASA
    tfe = transfer_free_energy(MJC, "urea"; model=MTRecord)
    tfe_sasa = transfer_free_energy(sasa_particles(RichardsUnitedAtomRadii, MJC), "urea"; model=MTRecord)
    @test tfe.tot ≈ tfe_sasa.tot

    # The paper has no equivalent extended-chain validation set for betaine (Table S3 is
    # urea-only), so these are regression checks against the model's own output, not
    # literature-validated targets.
    rec_m = MTRecordDenaturedModel(RN2)
    c_rec = mvalue(rec_m, "betaine").tot
    @test c_rec ≈ 0.9574 rtol = 0.05
    rec_m = MTRecordDenaturedModel(MJC)
    c_rec = mvalue(rec_m, "betaine").tot
    @test c_rec ≈ 0.4556 rtol = 0.05

    # Tests for arithmetic operations on TFEs
    tfe_n = transfer_free_energy(MJC, "urea"; model=MTRecord)
    tfe_test = 2 * tfe_n - tfe_n
    @test tfe_test.tot ≈ tfe_n.tot
    @test tfe_test.bb ≈ tfe_n.bb
    @test tfe_test.sc ≈ tfe_n.sc
    @test tfe_test.residue_contributions_bb ≈ tfe_n.residue_contributions_bb
    @test tfe_test.residue_contributions_sc ≈ tfe_n.residue_contributions_sc
    @test tfe_test.cosolvent == tfe_n.cosolvent
    tfe_n2 = transfer_free_energy(RN2, "urea"; model=MTRecord)
    @test_throws "number of residues" tfe_n2 + tfe_n
    tfe_n2 = transfer_free_energy(MJC, "betaine"; model=MTRecord)
    @test_throws "cosolvents" tfe_n2 + tfe_n

    # Test error path
    pdb = read_pdb(PDBTools.TESTPDB, "protein or resname TMAO")
    @test_throws "Creamer united atom" transfer_free_energy(pdb, "urea")
    @test_throws "residue TMAO" transfer_free_energy(pdb, "urea")

    # Test available cossolvents string
    using InteractiveUtils: subtypes
    models = replace.(string.(subtypes(PDBTools.MValueModel)), "PDBTools." => "")
    @test all(occursin(model_name,PDBTools._available_cosolvents()) for model_name in models)

    # Test modelname definitions
    @test PDBTools.modelname(MoeserHorinek) == "MoeserHorinek"
    @test PDBTools.modelname(AutonBolen) == "AutonBolen"
    @test PDBTools.modelname(MoeserHorinekApp) == "MoeserHorinekApp"
    @test PDBTools.modelname(Accessibility) == "Accessibility"
    @test PDBTools.modelname(MTRecord) == "MTRecord"

    # Model type returns
    @test PDBTools._model_type("MoeserHorinek") == MoeserHorinek
    @test PDBTools._model_type("AutonBolen") == AutonBolen
    @test PDBTools._model_type("MoeserHorinekApp") == MoeserHorinekApp
    @test PDBTools._model_type("Accessibility") == Accessibility
    @test PDBTools._model_type("MTRecord") == MTRecord
    @test_throws "Invalid MValueModel" PDBTools._model_type("ABC") 

    # Error if wrong atomic radii was provided
    s = sasa_particles(pdb)
    @test_throws "must use CreamerUnitedAtomRadii" transfer_free_energy(s, "urea")

end

@testitem "Record transfer model atom type mapping" begin
    using PDBTools
    using PDBTools: record_surface_type, model_combination_rule
    using ShowMethodTesting

    @test record_surface_type(Atom(resname="LEU", name="CB")) == :aliphatic_carbon
    @test record_surface_type(Atom(resname="PHE", name="CG")) == :aromatic_carbon
    @test record_surface_type(Atom(resname="SER", name="OG")) == :hydroxyl_oxygen
    @test record_surface_type(Atom(resname="ASN", name="OD1")) == :amide_oxygen
    @test record_surface_type(Atom(resname="ASP", name="OD2")) == :carboxylate_oxygen
    @test record_surface_type(Atom(resname="GLN", name="NE2")) == :amide_nitrogen
    @test record_surface_type(Atom(resname="HIS", name="ND1")) == :cationic_nitrogen
    @test record_surface_type(Atom(resname="TRP", name="NE1")) == :amide_nitrogen
    @test isnothing(record_surface_type(Atom(resname="MET", name="SD")))
    @test isnothing(record_surface_type(Atom(resname="CYS", name="SG")))
    @test_throws "unsupported element" record_surface_type(Atom(resname="ALA", name="P", pdb_element="P"))

    # Fallback for usuported types
    @test record_surface_type(Atom(resname="ARG", name="OAM")) == :amide_oxygen
    @test record_surface_type(Atom(resname="ARG", name="NAM")) == :amide_nitrogen

    # Registering the macro keywords a second time (e.g. on reload) must be a no-op.
    PDBTools._register_record_macro_keywords!()
    icn = findfirst(k -> k.name == "record_cationic_nitrogen", PDBTools.macro_keywords)
    popat!(PDBTools.macro_keywords, icn)
    PDBTools._register_record_macro_keywords!()
    @test findfirst(k -> k.name == "record_cationic_nitrogen", PDBTools.macro_keywords) !== nothing

    # Test macro-keyword integration with selection syntax.
    ats = [
        Atom(resname="ALA", name="N"),
        Atom(resname="ALA", name="CA"),
        Atom(resname="ALA", name="C"),
        Atom(resname="ALA", name="O"),
        Atom(resname="PHE", name="CG"),
        Atom(resname="SER", name="OG"),
        Atom(resname="ASN", name="OD1"),
        Atom(resname="ASP", name="OD2"),
        Atom(resname="GLN", name="NE2"),
        Atom(resname="HIS", name="ND1"),
    ]
    @test length(select(ats, "record_aliphatic_carbon and backbone")) == 2
    @test length(select(ats, "record_aromatic_carbon")) == 1
    @test length(select(ats, "record_hydroxyl_oxygen")) == 1
    @test length(select(ats, "record_amide_oxygen")) == 2
    @test length(select(ats, "record_carboxylate_oxygen")) == 1
    @test length(select(ats, "record_amide_nitrogen")) == 2
    @test length(select(ats, "record_cationic_nitrogen")) == 1

    # MTRecordDenaturedModel needs a real backbone (to build the extended chain from),
    # unlike the synthetic, disconnected `ats` used above.
    pep = read_pdb(PDBTools.TESTPDB, "protein and residue 1 to 5")
    m = MTRecordDenaturedModel(pep)
    @test parse_show(m) ≈ "MTRecordDenaturedModel wrapping a $(length(pep))-atom native/extended chain pair"
    @test length(m.native_chain) == length(pep)
    @test length(m.extended_chain) == length(pep)
    ram_ext = Ramachandran(m.extended_chain)
    @test all(x -> isapprox(x, 180.0; atol=1e-2) || isapprox(abs(x), 180.0; atol=1e-2), ram_ext.phi)
    @test all(x -> isapprox(x, 180.0; atol=1e-2) || isapprox(abs(x), 180.0; atol=1e-2), ram_ext.psi)

    @test isapprox(model_combination_rule(MTRecord, "urea", :aromatic_carbon), -0.5273f0; atol=1f-4)
    @test isapprox(model_combination_rule(MTRecord, "betaine", :amide_oxygen), 1.6590f0; atol=1f-4)

    # TetraEG and glycerol group interaction potentials are given directly in
    # cal mol⁻¹ molal⁻¹ Å⁻² in their source table, unlike the Guinn et al. cosolvents
    # above, so model_combination_rule must return them unscaled.
    @test isapprox(model_combination_rule(MTRecord, "glycerol", :hydroxyl_oxygen), 0.0305f0; atol=1f-4)
    @test isapprox(model_combination_rule(MTRecord, "glycerol", :cationic_nitrogen), -0.245f0; atol=1f-4)
    @test isapprox(model_combination_rule(MTRecord, "tetraeg", :carboxylate_oxygen), 3.59f0; atol=1f-4)
    @test isapprox(model_combination_rule(MTRecord, "tetraeg", :aliphatic_carbon), -0.349f0; atol=1f-4)

    tfe = transfer_free_energy(MTRecord, ats, "urea")
    @test tfe isa TransferFreeEnergy{MTRecord}
    @test tfe.nresidues == 7

    # Sanity check that the full pipeline runs for the new cosolvents.
    tfe_glycerol = transfer_free_energy(MTRecord, pep, "glycerol")
    @test tfe_glycerol isa TransferFreeEnergy{MTRecord}
    tfe_tetraeg = transfer_free_energy(MTRecord, pep, "tetraeg")
    @test tfe_tetraeg isa TransferFreeEnergy{MTRecord}
end

@testitem "tfe_dimer_tests" begin
    using PDBTools
    dir = @__DIR__
    dimer = read_pdb(joinpath(dir, "2rmm_dimer.pdb"))
    tfe = transfer_free_energy(dimer, "tmao"; model=Accessibility)
    @test tfe.tot ≈ 1.4553362

    # Remove chain identifiers
    dimer_nochain = copy.(dimer)
    for at in dimer_nochain
        at.chain = ""
    end
    tfe = transfer_free_energy(dimer_nochain, "tmao"; model=Accessibility)
    @test tfe.tot ≈ 1.4553362

    # Test single chains
    tfe1 = transfer_free_energy(select(dimer, "chain A"), "tmao"; model=Accessibility)
    tfe2 = transfer_free_energy(select(dimer_nochain, "residue <= 56"), "tmao"; model=Accessibility)
    @test tfe1.tot ≈ tfe2.tot
end
