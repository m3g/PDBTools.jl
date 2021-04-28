isacidic(atom::Atom) = 
  (i=which_natural_aminoacid(atom)) == 0 ? false : (natural_aminoacids[i].type == "Acidic")

isaliphatic(atom::Atom) = 
  (i=which_natural_aminoacid(atom)) == 0 ? false : (natural_aminoacids[i].type == "Aliphatic")

isaromatic(atom::Atom) = 
  (i=which_natural_aminoacid(atom)) == 0 ? false : (natural_aminoacids[i].type == "Aromatic")

isbasic(atom::Atom) =
  (i=which_natural_aminoacid(atom)) == 0 ? false : (natural_aminoacids[i].type == "Basic")

ischarged(atom::Atom) =
  (i=which_natural_aminoacid(atom)) == 0 ? false : (natural_aminoacids[i].charge != 0)

ishydrophobic(atom::Atom) =
  (i=which_natural_aminoacid(atom)) == 0 ? false : (natural_aminoacids[i].hydrophobic == true)

isneutral(atom::Atom) =
  (i=which_natural_aminoacid(atom)) == 0 ? false : (natural_aminoacids[i].charge == 0)

isnonpolar(atom::Atom) =
  (i=which_natural_aminoacid(atom)) == 0 ? false : (natural_aminoacids[i].polar == false)

ispolar(atom::Atom) =
  (i=which_natural_aminoacid(atom)) == 0 ? false : (natural_aminoacids[i].polar == true)

isprotein(atom::Atom; newres = nothing) = 
  (atom.resname == newres || which_natural_aminoacid(atom) != 0) ? true : false

const backbone_atoms = [ "N", "CA", "C", "O" ] 
isbackbone(atom::Atom) =
  which_natural_aminoacid(atom) == 0 ? false : (atom.name in backbone_atoms)

const not_side_chain_atoms = [ "N", "CA", "C", "O", "HN", "HA", "HT1", "HT2", "HT3" ]
issidechain(atom::Atom) =
  which_natural_aminoacid(atom) == 0 ? false : (! (atom.name in not_side_chain_atoms))

const water_residues = [ "HOH"  , "OH2"  ,
                         "TIP3" , "TIP3P",
                         "TIP4P", "TIP5P",
                         "TIP7P", "SPC"  ,
                         "SPCE" ] 
iswater(atom::Atom) = (findfirst(isequal(atom.resname),water_residues) != nothing) 

