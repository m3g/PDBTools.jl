#
# Check if structure is in mmCIF format and returns the list of 
# desired fields if that is the case
#

function check_mmCIF(file)

  local ifield, mmCIF, line

  mmCIF_fields = Indexes_mmCIF_fields()

  pdbfile = open(file,"r")

  mmCIF = false

  while ! eof(pdbfile)
    try 
      line = strip(readline(pdbfile))
      if length(line) < 5 
        continue
      end
    catch
      break
    end
    if line[1:5] == "loop_" 
      mmCIF = true
      ifield = 0
      while ! eof(pdbfile)
        try 
          line = strip(readline(pdbfile))
        catch
          break
        end
        if line[1:1] == "#" 
          continue
        end
        if ifield > 0 
          ifield = ifield + 1
          if occursin("_atom_site.label_atom_id",line) ; mmCIF_fields.name = ifield ; end 
          if occursin("_atom_site.id",line) ; mmCIF_fields.index = ifield ; end 
          if occursin("_atom_site.label_comp_id",line) ; mmCIF_fields.resname = ifield ; end 
          if occursin("_atom_site.label_asym_id",line) ; mmCIF_fields.chain = ifield ; end 
          if occursin("_atom_site.label_seq_id",line) ; mmCIF_fields.resnum = ifield ; end 
          if occursin("_atom_site.Cartn_x",line) ; mmCIF_fields.x = ifield ; end 
          if occursin("_atom_site.Cartn_y",line) ; mmCIF_fields.y = ifield ; end 
          if occursin("_atom_site.Cartn_z",line) ; mmCIF_fields.z = ifield ; end 
          if occursin("_atom_site.B_iso_or_equiv",line) ; mmCIF_fields.b = ifield ; end 
          if occursin("_atom_site.occupancy",line) ; mmCIF_fields.occup = ifield ; end 
          if occursin("ATOM",line) && occursin("HETATOM",line)
            break
          end
        end
        if occursin("_atom_site.group_PDB",line)
          ifield = 1
        end
      end
    end
  end
  close(pdbfile)

  if mmCIF
    if mmCIF_fields.name == 0 ||
       mmCIF_fields.resname == 0 ||
       mmCIF_fields.chain == 0 ||
       mmCIF_fields.resnum == 0 ||
       mmCIF_fields.x == 0 ||
       mmCIF_fields.y == 0 ||
       mmCIF_fields.z == 0 
      error(" ERROR: Could not find all necessary fields in mmCIF file. ")
    end
  end

  return mmCIF, mmCIF_fields

end



