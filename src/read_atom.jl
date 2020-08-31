#
# Function that reads atom information from PDB or mmCIF files
#

function read_atom(record :: String; 
                   mmCIF :: Bool = false, 
                   mmCIF_fields :: Indexes_mmCIF_fields = empty_struct(Indexes_mmCIF_fields))

  atom = MutableAtom()

  if length(record) < 6
    return nothing
  end

  if record[1:4] == "ATOM" || record[1:6] == "HETATM" 

    if ! mmCIF

      atom.name = strip(record[13:16])

      atom.resname = strip(record[17:21])
      resname = alternate_conformation(atom)
      if resname == nothing
        return nothing
      else
        atom.resname = resname
      end

      atom.chain = strip(record[22:22])
      if atom.chain == " "
        atom.chain = "0"
      end
      
      atom.index = 1
      try
        atom.index_pdb = parse_int(record[7:11])
      catch
        atom.index_pdb = 0
      end
      atom.resnum = parse_int(record[23:26])
      atom.x = parse(Float64,record[31:38])
      atom.y = parse(Float64,record[39:46])
      atom.z = parse(Float64,record[47:54])
      try
        atom.b = parse(Float64,record[61:66])
      catch
        atom.b = 0.
      end
      try 
        atom.occup = parse(Float64,record[56:60])
      catch
        atom.occup = 0.
      end
      atom.model = 1
      try 
        atom.segname = strip(record[73:76])
      catch
        atom.segname = ""
      end

    else # if mmCIF

      mmcif_data = split(record)
      atom.index = 1
      try
        atom.index_pdb = parse(Int64,mmcif_data[mmCIF_fields.index])
      catch
        atom.index_pdb = 0
      end
      atom.name = mmcif_data[mmCIF_fields.name]
      atom.resname = mmcif_data[mmCIF_fields.resname]
      resname = alternate_conformation(atom)
      if resname == nothing
        return nothing
      else
        atom.resname = resname
      end
      atom.chain = mmcif_data[mmCIF_fields.chain]
      try
        atom.resnum = parse(Int64,mmcif_data[mmCIF_fields.resnum])
      catch
        atom.resnum = 0
      end
      try
        atom.segname = mmcif_data[mmCIF_fields.segname]
      catch
        atom.segname = ""
      end
      atom.x = parse(Float64,mmcif_data[mmCIF_fields.x])
      atom.y = parse(Float64,mmcif_data[mmCIF_fields.y])
      atom.z = parse(Float64,mmcif_data[mmCIF_fields.z])
      atom.b = parse(Float64,mmcif_data[mmCIF_fields.b])
      atom.occup = parse(Float64,mmcif_data[mmCIF_fields.occup])
      atom.model = 1

    end

    return atom 

  else

    return nothing

  end
  
end
