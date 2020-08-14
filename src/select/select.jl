#
# Function to perform the most important selections that are used in solute-solvent
# analysis
#

# Main function: receives the atoms vector and a julia function to select

function select( atoms :: Union{Vector{PDBTools.Atom},Vector{PDBTools.MutableAtom}};
                 by=x->x==x)
  selected_atoms = typeof(atoms)(undef,0)
  for atom in atoms
    if by(atom)
      push!(selected_atoms,atom)
    end
  end
  return selected_atoms
end

# Given a selection string

function select( atoms :: Union{Vector{PDBTools.Atom},Vector{PDBTools.MutableAtom}}, 
                 selection :: String )
  query = parse_query(selection)
  return select(atoms, by = atom -> apply_query(query,atom))
end

#
# Return indexes only
#
function selindex( atoms :: Union{Vector{PDBTools.Atom},Vector{PDBTools.MutableAtom}};
                   by=x->x==x)
  indexes = Vector{Int64}(undef,0)
  for atom in atoms
    if by(atom)
      push!(indexes,atom.index)
    end
  end
  return indexes
end

function selindex( atoms :: Union{Vector{PDBTools.Atom},Vector{PDBTools.MutableAtom}}, 
                   selection :: String )
  query = parse_query(selection)
  return selindex(atoms, by = atom -> apply_query(query,atom))
end

# parse_query and apply_query are a very gentle contribution given by 
# CameronBieganek in https://discourse.julialang.org/t/parsing-selection-syntax/43632/9
# while explaining to me how to creat a syntex interpreter

function parse_query(selection)
  # disambiguate the "name" keyword
  s = replace(selection,"resname" => "RESNAME")
  try
    if occursin("or", s)
      (|, parse_query.(split(s, "or"))...)
    elseif occursin("and", s)
      (&, parse_query.(split(s, "and"))...)
    elseif occursin("not", s)
      rest = match(r".*not(.*)", s)[1]
      (!, parse_query(rest))

    # Index 
    elseif occursin("index =", s)
      k = parse(Int, match(r"index = ([0-9]*)", s)[1])
      a -> a.index == k
    elseif occursin("index < ", s)
      k = parse(Int, match(r"index < ([0-9]*)", s)[1])
      a -> a.index < k
    elseif occursin("index > ", s)
      k = parse(Int, match(r"index > ([0-9]*)", s)[1])
      a -> a.index > k
    elseif occursin("index <=", s)
      k = parse(Int, match(r"index <= ([0-9]*)", s)[1])
      a -> a.index <= k
    elseif occursin("index >=", s)
      k = parse(Int, match(r"index >= ([0-9]*)", s)[1])
      a -> a.index >= k

    # Resiue number
    elseif occursin("resnum =", s)
      k = parse(Int, match(r"resnum = ([0-9]*)", s)[1])
      a -> a.resnum == k
    elseif occursin("resnum < ", s)
      k = parse(Int, match(r"resnum < ([0-9]*)", s)[1])
      a -> a.resnum < k
    elseif occursin("resnum > ", s)
      k = parse(Int, match(r"resnum > ([0-9]*)", s)[1])
      a -> a.resnum > k
    elseif occursin("resnum <=", s)
      k = parse(Int, match(r"resnum <= ([0-9]*)", s)[1])
      a -> a.resnum <= k
    elseif occursin("resnum >=", s)
      k = parse(Int, match(r"resnum >= ([0-9]*)", s)[1])
      a -> a.resnum >= k

    # b factor
    elseif occursin("b =", s)
      b = parse(Float64, match(r"b = ([0-9]*)", s)[1])
      a -> a.b == k
    elseif occursin("b >", s)
      b = parse(Float64,match(r"b > ([.,0-9]*)", s)[1])
      a -> a.b > b
    elseif occursin("b <", s)
      b = parse(Float64,match(r"b < ([.,0-9]*)", s)[1])
      a -> a.b < b
    elseif occursin("b <=", s)
      b = parse(Float64,match(r"b <= ([.,0-9]*)", s)[1])
      a -> a.b <= b
    elseif occursin("b >=", s)
      b = parse(Float64,match(r"b >= ([.,0-9]*)", s)[1])
      a -> a.b >= b
    
    # occup
    elseif occursin("occup =", s)
      occup = parse(Float64, match(r"occup = ([0-9]*)", s)[1])
      a -> a.b == k
    elseif occursin("occup >", s)
      occup = parse(Float64,match(r"occup > ([.,0-9]*)", s)[1])
      a -> a.b > occup
    elseif occursin("occup <", s)
      occup = parse(Float64,match(r"occup < ([.,0-9]*)", s)[1])
      a -> a.b < occup
    elseif occursin("occup <=", s)
      occup = parse(Float64,match(r"occup <= ([.,0-9]*)", s)[1])
      a -> a.b <= occup 
    elseif occursin("occup >=", s)
      occup = parse(Float64,match(r"occup >= ([.,0-9]*)", s)[1])
      a -> a.b >= occup

    # Names etc.
    elseif occursin("name", s)
      name = match(r"name ([A-Z,0-9]*)", s)[1]
      a -> a.name == name
    elseif occursin("RESNAME", s)
      resname = match(r"RESNAME ([A-Z,0-9]*)", s)[1]
      a -> a.resname == resname
    elseif occursin("chain", s)
      chain = match(r"chain ([A-Z,0-9]*)", s)[1]
      a -> a.chain == chain
    elseif occursin("model", s)
      model = parse(Int,match(r"model ([0-9]*)", s)[1])
      a -> a.model == model
    elseif occursin("element", s)
      el = match(r"element ([A-Z]*)", s)[1]
      a -> element(a.name) == el

    # Special functions 
    elseif occursin("water", s)
      iswater
    elseif occursin("protein", s)
      isprotein
    elseif occursin("polar", s)
      ispolar
    elseif occursin("nonpolar", s)
      isnonpolar
    elseif occursin("basic", s)
      isbasic
    elseif occursin("acidic", s)
      isacidic
    elseif occursin("charged", s)
      ischarged
    elseif occursin("aliphatic", s)
      isaliphatic
    elseif occursin("aromatic", s)
      isaromatic
    elseif occursin("hydrophobic", s)
      ishydrophobic
    elseif occursin("neutral", s)
      isneutral
    elseif occursin("backbone", s)
      isbackbone
    elseif occursin("sidechain", s)
      issidechain

    # Select everything
    elseif occursin("all", s)
      a -> true

    # None found
    else
      parse_error()
    end

  # Error in the syntax
  catch err
    parse_error()
  end
end

function apply_query(q, a)
  if !(q isa Tuple)
    q(a)
  else
    f, args = Iterators.peel(q)
    f(apply_query.(args, Ref(a))...)
  end
end         

#
# Simple error message (https://discourse.julialang.org/t/a-julia-equivalent-to-rs-stop/36568/13)
#

struct NoBackTraceException
    exc::Exception
end

function Base.showerror(io::IO, ex::NoBackTraceException, bt; backtrace=true)
    Base.with_output_color(get(io, :color, false) ? Base.error_color() : :nothing, io) do io
        showerror(io, ex.exc)
    end
end

parse_error() = throw(NoBackTraceException(ErrorException("Error parsing selection. Use spaces, parenthesis not supported.")))

