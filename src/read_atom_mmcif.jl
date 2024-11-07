#=
data_all
loop_
_atom_site.group_PDB             # ATOM 
_atom_site.id                    # 1 index_pdb
_atom_site.type_symbol           # N symbol?  
_atom_site.label_atom_id         # N name
_atom_site.label_alt_id          # . 
_atom_site.label_comp_id         # VAL resname
_atom_site.label_asym_id         # A chain
_atom_site.label_entity_id       # 1 resnum
_atom_site.label_seq_id          # 1 ? 
_atom_site.pdbx_PDB_ins_code     # ? ? 
_atom_site.Cartn_x               # 6.204 x
_atom_site.Cartn_y               # 16.869 y
_atom_site.Cartn_z               # 4.854 z
_atom_site.occupancy             # 1.00 occup
_atom_site.B_iso_or_equiv        # 49.05 beta
_atom_site.pdbx_formal_charge    # ? charge
_atom_site.auth_seq_id           # 1
_atom_site.auth_comp_id          # VAL
_atom_site.auth_asym_id          # A
_atom_site.auth_atom_id          # N
_atom_site.pdbx_PDB_model_num    # 1 model
ATOM   1    N  N   . VAL A 1 1   ? 6.204   16.869  4.854   1.00 49.05 ? 1   VAL A N   1
ATOM   2    C  CA  . VAL A 1 1   ? 6.913   17.759  4.607   1.00 43.14 ? 1   VAL A CA  1
ATOM   3    C  C   . VAL A 1 1   ? 8.504   17.378  4.797   1.00 24.80 ? 1   VAL A C   1
ATOM   4    O  O   . VAL A 1 1   ? 8.805   17.011  5.943   1.00 37.68 ? 1   VAL A O   1
ATOM   5    C  CB  . VAL A 1 1   ? 6.369   19.044  5.810   1.00 72.12 ? 1   VAL A CB  1
ATOM   6    C  CG1 . VAL A 1 1   ? 7.009   20.127  5.418   1.00 61.79 ? 1   VAL A CG1 1
ATOM   7    C  CG2 . VAL A 1 1   ? 5.246   18.533  5.681   1.00 80.12 ? 1   VAL A CG2 1 
=#

# read atom from mmCIF file
function read_atom_mmCIF(record::String, mmCIF_fields::indices_mmCIF_fields=indices_mmCIF_fields())
    if !startswith(record, r"ATOM|HETATM")
        return nothing
    end
    atom = Atom()
    mmcif_data = split(record)
    atom.index = 1
    atom.index_pdb = parse_number(Int, mmcif_data[mmCIF_fields.index])
    atom.name = mmcif_data[mmCIF_fields.type_symbol]
    atom.resname = mmcif_data[mmCIF_fields.label_comp_id]
    atom.chain = mmcif_data[mmCIF_fields.chain]
    atom.resnum = parse_number(Int, mmcif_data[mmCIF_fields.label_asym_id])
    atom.segname = mmcif_data[mmCIF_fields.segname]
    try
        atom.segname = mmcif_data[mmCIF_fields.segname]
    catch
        atom.segname = ""
    end
    atom.x = parse(Float64, mmcif_data[mmCIF_fields.x])
    atom.y = parse(Float64, mmcif_data[mmCIF_fields.y])
    atom.z = parse(Float64, mmcif_data[mmCIF_fields.z])
    atom.beta = parse(Float64, mmcif_data[mmCIF_fields.beta])
    atom.occup = parse(Float64, mmcif_data[mmCIF_fields.occup])
    atom.model = 1
    return atom
end