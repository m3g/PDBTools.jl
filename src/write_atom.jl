#
# Function that returns an ATOM line in PDB format
#

# requires Printf

function write_atom(atom :: MutableAtom)
  imutt_atom = Atom(atom)
  return write_atom(imutt_atom)
end

function write_atom(atom :: Atom)

  #ATOM      2  CA  GLY A   1      -1.774   6.778  32.054  1.00  0.08           C
  #COLUMNS        DATA  TYPE    FIELD        DEFINITION
  #-------------------------------------------------------------------------------------
  # 1 -  6        Record name   "ATOM  "
  # 7 - 11        Integer       serial       Atom  serial number.
  #13 - 16        Atom          name         Atom name.
  #17             Character     altLoc       Alternate location indicator.
  #18 - 20        Residue name  resName      Residue name.
  #22             Character     chainID      Chain identifier.
  #23 - 26        Integer       resSeq       Residue sequence number.
  #27             AChar         iCode        Code for insertion of residues.
  #31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
  #39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
  #47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
  #55 - 60        Real(6.2)     occupancy    Occupancy.
  #61 - 66        Real(6.2)     tempFactor   Temperature  factor.
  #77 - 78        LString(2)    element      Element symbol, right-justified.
  #79 - 80        LString(2)    charge       Charge  on the atom.

  name = strip(atom.name)
  l = length(name)
  if l == 1
    name = "  $(atom.name) "
  elseif l == 2
    name = " $(atom.name) "
  elseif l == 3
    name = " $(atom.name)"
  else
    name = atom.name
  end
  if isprotein(atom)
    resname = strip(atom.resname)
    l = length(resname)
    if l == 3
      resname = " $(atom.resname) "
    elseif l == 4
      resname = "$(atom.resname) "
    end
  else
    resname = strip(atom.resname)
    l = length(resname)
    if l == 1
      resname = " $(atom.resname)   "
    elseif l == 2
      resname = " $(atom.resname)  "
    elseif l == 3
      resname = " $(atom.resname) "
    elseif l == 4
      resname = " $(atom.resname)"
    else
      resname = atom.resname
    end
  end

  if abs(atom.x) > 999.999 || 
     abs(atom.y) > 999.999 || 
     abs(atom.y) > 999.999  
    println("Warning: coordinates of atom $atom.index"*
            "(atom.name atom.resname$atom.resnum atom.chain)"*
            "greater than 999.999 and do not fit in PDB fixed format.")
  end
   
  if atom.index < 100000
    line = @sprintf("%-6s%5i%1s%4s%4s%1s%4i%4s%8.3f%8.3f%8.3f%6.2f%6.2f",
                     "ATOM",
                     atom.index," ",name,resname,atom.chain,atom.resnum,"    ",
                     atom.x,atom.y,atom.z,atom.occup,atom.b)
  else # Prints hexadecimal code for atom index
    line = @sprintf("%-6s%5x%1s%4s%4s%1s%4i%4s%8.3f%8.3f%8.3f%6.2f%6.2f",
                     "ATOM",
                     atom.index," ",name,resname,atom.chain,atom.resnum,"    ",
                     atom.x,atom.y,atom.z,atom.occup,atom.b)
  end

  return line
end


