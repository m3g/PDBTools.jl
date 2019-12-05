#
# Function that returns an ATOM line in PDB format
#

# requires Printf

function write_atom(atom :: ReadAtom)
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

  l = length(atom.name)
  if l == 1
    name = "  $(atom.name) "
  elseif l == 2
    name = " $(atom.name) "
  elseif l == 3
    name = " $(atom.name)"
  else
    name = atom.name
  end
  line = @sprintf("%-6s%5i%1s%4s%1s%-3s%1s%1s%4i%4s%8.3f%8.3f%8.3f%6.2f%6.2f",
                   "ATOM",
                   atom.index," ",name," ",atom.resname," ",atom.chain,atom.resnum,"    ",
                   atom.x,atom.y,atom.z,atom.occup,atom.b)

  return line
end


