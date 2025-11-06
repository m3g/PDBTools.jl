@testitem "mvalue" begin
    using PDBTools
    using ShowMethodTesting

    dir = @__DIR__
    MJC_native = read_pdb(joinpath(dir, "1MJC_native.pdb"), "protein")
    MJC_desnat = read_pdb(joinpath(dir, "1MJC_straight.pdb"), "protein")

    #
    # AutonBolen model
    #
    n_sasa = sasa_particles(MJC_native)
    d_sasa = sasa_particles(MJC_desnat)
    r_1MJC = mvalue(n_sasa, d_sasa, "urea")
    @test isapprox(r_1MJC.tot, -0.786; atol=1e-2)
    @test isapprox(r_1MJC.bb, -0.852; atol=1e-2)
    @test isapprox(r_1MJC.sc, 0.066; atol=1e-2)
    # Confirm that now cosolvent selection is not case-sensitive
    r_1MJC = mvalue(n_sasa, d_sasa, "Urea")
    @test isapprox(r_1MJC.tot, -0.786; atol=1e-2)

    # Test show method
    @test parse_show(r_1MJC; repl=Dict(r"PDBTools." => "")) ≈ """
        MValue{AutonBolen} - 69 residues.
            Total m-value: -0.7856683 kcal mol⁻¹
            Backbone contributions: -0.8520461 kcal mol⁻¹
            Side-chain contributions: 0.06637777 kcal mol⁻¹
        """

    # Without H
    n_sasa = sasa_particles(select(MJC_native, "protein and not element H"))
    d_sasa = sasa_particles(select(MJC_desnat, "protein and not element H"))
    r_1MJC = mvalue(n_sasa, d_sasa, "urea")
    @test isapprox(r_1MJC.tot, -1.546; atol=1e-2)
    @test isapprox(r_1MJC.bb, -1.600; atol=1e-2)
    @test isapprox(r_1MJC.sc, 0.0539; atol=1e-2)

    r_1MJC = mvalue(n_sasa, d_sasa, "tmao")
    @test isapprox(r_1MJC.tot, 3.073; atol=1e-2)
    @test isapprox(r_1MJC.bb, 3.692; atol=1e-2)
    @test isapprox(r_1MJC.sc, -0.619; atol=1e-2)

    #
    # MoeserHorinek model
    # 
    n_sasa = sasa_particles(MJC_native)
    d_sasa = sasa_particles(MJC_desnat)
    r_1MJC = mvalue(n_sasa, d_sasa, "urea"; model=MoeserHorinek)
    @test isapprox(r_1MJC.tot, -0.937; atol=1e-2)
    @test isapprox(r_1MJC.bb, -0.391; atol=1e-2)
    @test isapprox(r_1MJC.sc, -0.546; atol=1e-2)

    # Without H
    n_sasa = sasa_particles(select(MJC_native, "protein and not element H"))
    d_sasa = sasa_particles(select(MJC_desnat, "protein and not element H"))
    r_1MJC = mvalue(n_sasa, d_sasa, "urea"; model=MoeserHorinek)
    @test isapprox(r_1MJC.tot, -1.222; atol=1e-2)
    @test isapprox(r_1MJC.bb, -0.757; atol=1e-2)
    @test isapprox(r_1MJC.sc, -0.465; atol=1e-2)

    # Same definitions as Gromacs for side chain and backbone
    n_sasa = sasa_particles(select(MJC_native, "protein and not element H"))
    d_sasa = sasa_particles(select(MJC_desnat, "protein and not element H"))
    r_1MJC = mvalue(n_sasa, d_sasa, "urea";
        model=MoeserHorinek,
        backbone=at -> name(at) in ("N", "CA", "C", "O"),
        sidechain=at -> !(name(at) in ("N", "CA", "C", "O")),
    )
    @test isapprox(r_1MJC.tot, -1.227; atol=1e-2)
    @test isapprox(r_1MJC.bb, -0.757; atol=1e-2)
    @test isapprox(r_1MJC.sc, -0.470; atol=1e-2)

    # Provide a selection
    r_1MJC = mvalue(n_sasa, d_sasa, "urea"; sel="acidic")
    @test isapprox(r_1MJC.tot, -0.037; atol=1e-2)
    @test isapprox(r_1MJC.bb, -0.104; atol=1e-2)
    @test isapprox(r_1MJC.sc, 0.067; atol=1e-2)

    # Input file with non-protein atoms
    pdb = read_pdb(PDBTools.TESTPDB, "protein or resname TMAO")
    s1 = sasa_particles(pdb)
    s2 = sasa_particles(select(pdb, "protein"))
    # This should work:
    m = mvalue(s1, s2, "urea"; sel="protein")
    @test m.tot ≈ 0.0 atol = 1e-2
    m = mvalue(s2, s1, "urea"; sel="protein")
    @test m.tot ≈ 0.0 atol = 1e-2
    # But these should error:
    @test_throws "same number of residues" mvalue(s2, s1, "urea") # different number of residues
    s2 = sasa_particles(pdb)
    @test_throws "non-protein residue" mvalue(s2, s1, "urea") # different number of residues
    s1 = sasa_particles(select(pdb, "protein"))
    pdb2 = copy.(pdb)
    rs = collect(eachresidue(pdb2))
    for at in rs[2]
        at.resname = "ASP"
    end
    s2 = sasa_particles(select(pdb2, "protein"))
    @test_throws "same type" mvalue(s2, s1, "urea") # different number of residues

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
    @test isapprox(r_1MJC.tot, -1.02; rtol=1e-1)
    @test isapprox(r_1MJC.bb, -0.389; rtol=1e-1)
    @test isapprox(r_1MJC.sc, -0.627; rtol=1e-1)

    if !isnothing(gmx)
        sasa_1MJC = gmx_delta_sasa_per_restype(;
            native_pdb=joinpath(dir, "1MJC_native.pdb"),
            desnat_pdb=joinpath(dir, "1MJC_straight.pdb"),
        )
        r_1MJC = mvalue_delta_sasa(; atoms=MJC_native, sasas=sasa_1MJC)
        @test isapprox(r_1MJC.tot, -1.24; rtol=1e-1)
        @test isapprox(r_1MJC.bb, -0.777; rtol=1e-1)
        @test isapprox(r_1MJC.sc, -0.463; rtol=1e-1)
        sasa_1MJC = gmx_delta_sasa_per_restype(;
            native_pdb=joinpath(dir, "1MJC_native.pdb"),
            desnat_pdb=joinpath(dir, "1MJC_straight.pdb"),
            ignore_hydrogen=false,
        )
        r_1MJC = mvalue_delta_sasa(; atoms=MJC_native, sasas=sasa_1MJC)
        @test isapprox(r_1MJC.tot, -1.01; rtol=1e-1)
        @test isapprox(r_1MJC.bb, -0.342; rtol=1e-1)
        @test isapprox(r_1MJC.sc, -0.668; rtol=1e-1)
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