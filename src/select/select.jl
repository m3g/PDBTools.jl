#
# Function to perform the most important selections that are used in solute-solvent
# analysis
#

function select( atoms :: Vector{PDBTools.Atom}, selection :: String )
  query = parse_query(selection)
  selection = Vector{PDBTools.Atom}(undef,0)
  for atom in atoms
    if apply_query(query,atom) 
      push!(selection,atom)
    end
  end
  return selection
end

function selindex( atoms :: Vector{PDBTools.Atom}, selection :: String )
  query = parse_query(selection)
  selection = Vector{Int64}(undef,0)
  for atom in atoms
    if apply_query(query,atom) 
      push!(selection,atom.index)
    end
  end
  return selection
end

# parse_query and apply_query are a very gentle contribution given by 
# CameronBieganek in https://discourse.julialang.org/t/parsing-selection-syntax/43632/9
# while explaining to me how to creat a syntex interpreter

function parse_query(s)
  try
    if occursin("or", s)
      (|, parse_query.(split(s, "or"))...)
    elseif occursin("and", s)
      (&, parse_query.(split(s, "and"))...)
    elseif occursin("not", s)
      rest = match(r".*not(.*)", s)[1]
      (!, parse_query(rest))
    elseif occursin("index < ", s)
      k = parse(Int, match(r"index < ([0-9]*)", s)[1])
      a -> a.index < k
    elseif occursin("index > ", s)
      k = parse(Int, match(r"index > ([0-9]*)", s)[1])
      a -> a.index > k
    elseif occursin("index =", s)
      k = parse(Int, match(r"index = ([0-9]*)", s)[1])
      a -> a.index == k
    elseif occursin("index >=", s)
      k = parse(Int, match(r"index >= ([0-9]*)", s)[1])
      a -> a.index >= k
    elseif occursin("index <=", s)
      k = parse(Int, match(r"index <= ([0-9]*)", s)[1])
      a -> a.index <= k
    elseif occursin("resname", s)
      resname = match(r"resname ([A-Z,0-9]*)", s)[1]
      a -> a.resname == resname
    elseif occursin("resnum", s)
      resnum = parse(Int, match(r"resnum ([0-9]*)", s)[1])
      a -> a.resnum == resname
    elseif occursin("residue", s)
      residue = parse(Int, match(r"residue ([0-9]*)", s)[1])
      a -> a.residue == residue
    elseif occursin("water", s)
      iswater
    elseif occursin("protein", s)
      isprotein
    else
      parse_error()
    end
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

