#= 

Isolated ASA values are from the Supporting Table 2 of https://doi.org/10.1073/pnas.0507053102
(https://www.pnas.org/doi/suppl/10.1073/pnas.0507053102/suppl_file/07053table2.pdf)

=#
const isolated_ASA = Dict{String,Tuple{Float32,Float32}}( 
                # BB      SC   (Å^2)
    "ALA"	=> (46.2,	71.9),
    "PHE"	=> (38.4,	184.4),
    "LEU"	=> (35.3,	157.8),
    "ILE"	=> (30.9,	150.1),
    "VAL"	=> (36.1,	128.4),
    "PRO"	=> (35.6,	111.0),
    "MET"	=> (38.6,	164.8),
    "TRP"	=> (37.4,	228.9),
    "GLY"	=> (88.1,	0),
    "SER"	=> (44.0,	85.8),
    "THR"	=> (37.9,	114.6),
    "TYR"	=> (38.7,	198.1),
    "GLN"	=> (37.8,	155.4),
    "ASN"	=> (40.2,	125.3),
    "ASP"	=> (40.5,	118.2),
    "GLU"	=> (37.8,	148.4),
    "HIS"	=> (40.4,	162.1),
    "LYS"	=> (38.7,	187.1),
    "ARG"	=> (39.1,	216.9),
    "CYS"	=> (42.6,	103.5),
) 