#
# Function to perform the most important selections that are used in solute-solvent
# analysis
#

# Main function: receives the atoms vector and a julia function to select

function select( atoms :: Vector{Atom}; by=all)
  selected_atoms = Atom[]
  for atom in atoms
    if by(atom)
      push!(selected_atoms,atom)
    end
  end
  return selected_atoms
end

# Given a selection string

function select( atoms :: Vector{Atom}, selection :: String )
  query = parse_query(selection)
  return select(atoms, by = atom -> apply_query(query,atom))
end

#
# Return indexes only
#
function selindex( atoms :: Vector{Atom}; by=all)
  indexes = Vector{Int64}(undef,0)
  for atom in atoms
    if by(atom)
      push!(indexes,atom.index)
    end
  end
  return indexes
end

function selindex( atoms :: Vector{Atom}, selection :: String )
  query = parse_query(selection)
  return selindex(atoms, by = atom -> apply_query(query,atom))
end

# Comparison operators

const operators = ( " = " =>  (x,y) -> isequal(x,y),
                    " < " =>  (x,y) -> isless(x,y),
                    " > " =>  (x,y) -> isless(y,x), 
                    " <= " => (x,y) -> (! isless(y,x)),
                    " >= " => (x,y) -> (! isless(x,y)) ) 

using Parameters

#
# Numerical syntax keywords
#

struct NumericalKeyword
  name :: String
  symbol :: Symbol
  operations :: Tuple
end

function parse_keyword(key::NumericalKeyword,str,op) 
  return parse(Int, match(key.name*op*r"([0-9]*)", str)[1])   
end

#
# String syntax keywords
#

struct StringKeyword
  name :: String
  symbol :: Symbol
  operations :: Tuple
end

function parse_keyword(key::StringKeyword,str,op) 
  return match(key.name*op*r"([A-Z,0-9]*)", str)[1] 
end

#
# Function that returns the anonymous functions that evaluates if a condition
# is satisfied or not for a given atom, given numerical or string keywords
#

function getfunc(key::T, str) where T <: Union{NumericalKeyword,StringKeyword}
  @unpack name, operations, symbol = key
  for op in operations
    if occursin(name*op.first, str)
      val = parse_keyword(key,str,op.first)
      return atom -> op.second(getfield(atom,symbol),val)
    end
  end
  # If got here, no operator was found, then we assume is is an implicit equal
  val = parse_keyword(key,str," ")
  return atom -> isequal(getfield(atom,symbol),val)
end

numerical_keywords = [ NumericalKeyword("index", :index, operators), 
                       NumericalKeyword("index_pdb",  :index_pdb, operators),
                       NumericalKeyword("resnum", :resnum, operators),
                       NumericalKeyword("residue", :residue, operators),
                       NumericalKeyword("b", :b, operators),
                       NumericalKeyword("occup", :occup, operators),
                       NumericalKeyword("model", :model, operators),
                     ]

string_keywords = [ StringKeyword("name", :name, operators), 
                    StringKeyword("segname", :segname, operators), 
                    StringKeyword("resname", :resname, operators), 
                    StringKeyword("chain", :chain, operators), 
                    StringKeyword("element", :element, operators), 
                  ]

#
# Special functions keywords
#

struct MacroKeyword
  name :: String
  fname :: Function
end
getfunc(key::MacroKeyword,str) = key.fname

macro_keywords = [ MacroKeyword("water", iswater), 
                   MacroKeyword("protein", isprotein), 
                   MacroKeyword("polar", ispolar), 
                   MacroKeyword("nonpolar", isnonpolar), 
                   MacroKeyword("basic", isbasic), 
                   MacroKeyword("acidic", isacidic), 
                   MacroKeyword("charged", ischarged), 
                   MacroKeyword("aliphatic", isaliphatic), 
                   MacroKeyword("aromatic", isaromatic), 
                   MacroKeyword("hydrophobic", ishydrophobic), 
                   MacroKeyword("neutral", isneutral), 
                   MacroKeyword("backbone", isbackbone), 
                   MacroKeyword("sidechain", issidechain), 
                   MacroKeyword("all", isequal), 
                 ]

#
# Remove trailing and multiple spaces
#

function remove_spaces(str)
  str = split(str)
  s = str[1]
  for i in 2:length(str)
    s = s*' '*str[i]
  end
  return " "*s*" "
end

# parse_query and apply_query are a very gentle contribution given by 
# CameronBieganek in https://discourse.julialang.org/t/parsing-selection-syntax/43632/9
# while explaining to me how to creat a syntex interpreter

function parse_query(selection)
  # Remove spaces
  s = remove_spaces(selection) 

  # Parse syntax
  try
    if occursin(" or ", s)
      (|, parse_query.(split(s, "or"))...)
    elseif occursin(" and ", s)
      (&, parse_query.(split(s, "and"))...)
    elseif occursin(" not ", s)
      rest = match(r".*not(.*)", s)[1]
      (!, parse_query(rest))

    # keywords 
    else
      for key in numerical_keywords
        occursin(" "*key.name*" ", s) && return getfunc(key,s)
      end
      for key in string_keywords
        occursin(" "*key.name*" ", s) && return getfunc(key,s)
      end
      for key in macro_keywords
        occursin(" "*key.name*" ", s) && return getfunc(key,s)
      end
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

