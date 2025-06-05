#
# Function to perform the most important selections that are used in solute-solvent
# analysis
#

using Parameters
export select, selindex
export Select, @sel_str

"""
    Select

This structure acts a function when used within typical julia filtering functions, 
by converting a string selection into a call to query call. 

# Example

Using a string to select the CA atoms of the first residue:

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.TESTPDB, "protein");

julia> findfirst(Select("name CA"), atoms)
5

julia> filter(Select("name CA and residue 1"), atoms)
   Vector{Atom{Nothing}} with 1 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       5   CA     ALA     A        1        1   -8.483  -14.912   -6.726  1.00  0.00     1    PROT         5

```

"""
struct Select{Q} <: Function
    query_string::String
    query::Q
end
function Select(query_string::AbstractString)
    query = parse_query(query_string)
    return Select(query_string, query)
end
(s::Select)(at) = apply_query(s.query, at)
Base.show(io::IO, ::MIME"text/plain", s::Select) = print(io, """Select("$(s.query_string)")""")

#
# Parse selection string allowing interpolation in sel macro:
# https://discourse.julialang.org/t/str-string-interpolation/125766/11?u=lmiq
#
_select(args...) = Select(string(args...))
macro sel_str(s)
    ex = Expr(:call, GlobalRef(PDBTools, :_select))
    i = firstindex(s)
    buf = IOBuffer(maxsize=ncodeunits(s))
    while i <= ncodeunits(s)
        c = @inbounds s[i]
        i = nextind(s, i)
        if c === '$'
            position(buf) > 0 && push!(ex.args, String(take!(buf)))
            val, i = Meta.parse(s, i; greedy=false)
            Meta.isexpr(val, :incomplete) && error(val.args[1])
            val !== nothing && push!(ex.args, val)
        else
            print(buf, c)
        end
    end
    position(buf) > 0 && push!(ex.args, String(take!(buf)))
    return esc(ex)
end

# Function that returns true for all atoms: the default selection
all(atoms) = true

@testitem "Select and sel_str show" begin
    using PDBTools
    atoms = read_pdb(PDBTools.TESTPDB, "protein")
    sel = Select("name CA and residue 1")
    buff = IOBuffer()
    show(buff, MIME"text/plain"(), sel)
    @test String(take!(buff)) == """Select("name CA and residue 1")"""
    show(buff, MIME"text/plain"(), sel"name CA and residue 1")
    @test String(take!(buff)) == """Select("name CA and residue 1")"""
end

"""
    select(atoms::AbstractVector{<:Atom}, by::String)

Selects atoms from a vector of atoms using a string query, or a function.

"""
function select end

# Main function: receives the atoms vector and a julia function to select
select(set::AbstractVector{<:Atom}, by::String) = filter(Select(by), set)
select(set::AbstractVector{<:Atom}, by::Function) = filter(by, set)

# Select indices of atoms (this is identical to findall)
selindex(set::AbstractVector{<:Atom}, by::String) = findall(Select(by), set)
selindex(set::AbstractVector{<:Atom}, by::Function) = findall(by, set)

# These two methods probably will be deprecated
function select(set::AbstractVector; by=all)
    @warn begin
        """\n
        The `select` function with the keyword argument `by` will be deprecated. 
        Use `select(atoms, function)` instead. Or simply `filter(function, atoms)`.

        """
    end _file = nothing _line = nothing
    filter(by, set)
end
function selindex(set::AbstractVector; by=all)
    @warn begin
        """\n
        The `selindex` function with the keyword argument `by` will be deprecated. 
        Use `findall(function, atoms)` instead."
            
        """
    end _file = nothing _line = nothing
    findall(by, set)
end

# Comparison operators

const operators = (
    "=" => (x, y) -> isequal(x, y),
    "<" => (x, y) -> isless(x, y),
    ">" => (x, y) -> isless(y, x),
    "<=" => (x, y) -> (!isless(y, x)),
    ">=" => (x, y) -> (!isless(x, y)),
)

#
# Keywords
#

struct Keyword
    ValueType::Type
    name::String
    field::Symbol
    operators::Tuple
end

function (key::Keyword)(s::AbstractVector{<:AbstractString})
    @unpack name, field, operators = key
    for op in operators
        if (i = has_key(op.first, s)) > 0
            return el -> op.second(getfield(el, field), parse_to_type(key, s[i+1]))
        end
    end
    # If no operator was found, assume that `=` was intended
    i = has_key(name, s)
    return el -> isequal(getfield(el, field), parse_to_type(key, s[i+1]))
end

#=
    FunctionalKeyword{T}

This is a structure that will store a keyword that depends on an external function
requiring an operator and an argument. 

## Example:

```
element_keyword = FunctionalKeyword(String,      
                                    "element",
                                    element,
                                    ("=",isequal))

```
will define a keyword "element" to be used as `element C`, which will return
`true` if there is an `element` function such that `element(atom) == C`.

=#
struct FunctionalKeyword{FunctionType}
    ValueType::Type
    name::String
    by::FunctionType
    operators::Tuple
end

function (key::FunctionalKeyword)(s::AbstractVector{<:AbstractString})
    @unpack name, by, operators = key
    for op in operators
        if (i = has_key(op.first, s)) > 0
            return el -> op.second(by(el), parse_to_type(key, s[i+1]))
        end
    end
    # If no operator was found, assume that `=` was intended
    i = has_key(name, s)
    return el -> isequal(by(el), parse_to_type(key, s[i+1]))
end

#
# Macro keywords (functions without parameters)
#
struct MacroKeyword{FunctionType}
    name::String
    by::FunctionType
end

function (key::MacroKeyword)(s::AbstractVector{<:AbstractString})
    return key.by
end

#=
    parse_to_type(key::Keyword, val::String)

Tries to parse `val` into the type of value expected by `key.ValueType`. 

=#
function parse_to_type(key::Union{Keyword,FunctionalKeyword}, val)
    if key.ValueType == String
        return val
    end
    try
        val = parse(key.ValueType, val)
        return val
    catch
        parse_error(
            """\n
                Could not parse $val for keyword $(key.name), expected $(key.ValueType)

            """,
        )
    end
end

#
# Keywords for PDBTools
#
keywords = [
    Keyword(Int, "index", :index, operators),
    Keyword(Int, "index_pdb", :index_pdb, operators),
    Keyword(Int, "resnum", :resnum, operators),
    Keyword(Int, "residue", :residue, operators),
    Keyword(Float64, "beta", :beta, operators),
    Keyword(Float64, "occup", :occup, operators),
    Keyword(Int, "model", :model, operators),
    Keyword(String, "name", :name, operators),
    Keyword(String, "segname", :segname, operators),
    Keyword(String, "resname", :resname, operators),
    Keyword(String, "chain", :chain, operators),
]

macro_keywords = [
    MacroKeyword("water", iswater),
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
    MacroKeyword("all", all),
]

functional_keywords = [
    FunctionalKeyword(String, "element", element, operators)
]

#
# parse_query and apply_query are a very gentle contribution given by 
# CameronBieganek in https://discourse.julialang.org/t/parsing-selection-syntax/43632/9
# while explaining to me how to creat a syntex interpreter
#
#=
    has_key(key::String, s::AbstractVector{<:AbstractString})

Returns the first index of the vector `s` in which where `key` is found, or 0. 

## Example:

```julia-repl

julia> PDBTools.has_key("or",["name","CA","or","index","1"])
3

julia> PDBTools.has_key("and",["name","CA","or","index","1"])
0

```

=#
function has_key(key::String, s::AbstractVector{<:AbstractString})
    i = findfirst(isequal(key), s)
    if isnothing(i)
        0
    else
        i
    end
end

#=
    parse_query(selection:String)

Calls `parse_query_vector` after splitting the selection string.

=#
function parse_query(selection::String) 
    s = replace(selection, "(" => " ( ", ")" => " ) ")
    return parse_query_vector(split(s))
end

function apply_query(q, a)
    if !(q isa Tuple)
        q(a)
    else
        f, args = Iterators.peel(q)
        f(apply_query.(args, Ref(a))...)
    end
end

parse_error(str) = throw(ArgumentError(str))

#
# Obs: the following code were generated by Gemini 2.5-Pro, with modifications, 
# and then tested. 
#

# New helper functions
function is_operator(token::AbstractString)
    return token == "and" || token == "or" || token == "not"
end

function is_fully_enclosed(tokens::AbstractVector{<:AbstractString})
    if length(tokens) < 2 || !(tokens[begin] == "(" && tokens[end] == ")")
        return false
    end
    level = 0
    # Check if the first '(' matches the last ')' without level becoming zero in between
    # for any token except the last one.
    for i in firstindex(tokens):(lastindex(tokens)-1)
        if tokens[i] == "("
            level += 1
        elseif tokens[i] == ")"
            level -= 1
            if level == 0 # Closed too early, means not fully enclosed by the outermost pair
                 return false
            end
        end
    end
    # After iterating up to tokens[end-1], level should be 1 if tokens[begin] was '('
    # and it correctly matches tokens[end]. If level is not 1, it means mismatched parentheses within.
    return level == 1
end

function find_operator_at_level_zero(op_str::String, tokens::AbstractVector{<:AbstractString})
    level = 0
    # Find first occurrence from left to right (maintaining current style)
    for i in eachindex(tokens)
        if tokens[i] == "("
            level += 1
        elseif tokens[i] == ")"
            level -= 1
            if level < 0
                parse_error("Mismatched parentheses: too many closing parentheses.")
            end
        elseif tokens[i] == op_str && level == 0
            return i
        end
    end
    if level != 0
        parse_error("Mismatched parentheses: not enough closing parentheses.")
    end
    return 0 # Not found at level zero
end

# Modified parse_query_vector
function parse_query_vector(s_vec_const::AbstractVector{<:AbstractString})
    s_vec = s_vec_const # Operate on slices or copies, not modifying original array passed around

    if isempty(s_vec)
        parse_error("Empty query segment.")
    end

    # Handle expressions fully enclosed in matching parentheses
    # e.g. "(A and B)" should be parsed by parsing "A and B"
    temp_s_vec = s_vec # Use a temporary variable for iterative stripping
    while length(temp_s_vec) > 1 && temp_s_vec[begin] == "(" && temp_s_vec[end] == ")" && is_fully_enclosed(temp_s_vec)
        temp_s_vec = temp_s_vec[begin+1:end-1]
        if isempty(temp_s_vec)
            parse_error("Empty parentheses in query: '()'")
        end
    end
    s_vec = temp_s_vec # Assign the stripped version back

    # Operator precedence: OR, then AND, then NOT (as in original code for splitting)
    # Find 'or' not within parentheses
    if (i = find_operator_at_level_zero("or", s_vec)) > 0
        left_tokens = s_vec[begin:i-1]
        right_tokens = s_vec[i+1:end]
        if isempty(left_tokens) || isempty(right_tokens)
            parse_error("Syntax error near 'or'. Missing operand.")
        end
        return (|, parse_query_vector(left_tokens), parse_query_vector(right_tokens))

    elseif (i = find_operator_at_level_zero("and", s_vec)) > 0
        left_tokens = s_vec[begin:i-1]
        right_tokens = s_vec[i+1:end]
        if isempty(left_tokens) || isempty(right_tokens)
            parse_error("Syntax error near 'and'. Missing operand.")
        end
        return (&, parse_query_vector(left_tokens), parse_query_vector(right_tokens))

    elseif s_vec[begin] == "not"
        if length(s_vec) == 1
            parse_error("Syntax error near 'not'. Missing operand.")
        end
        remaining_tokens = s_vec[begin+1:end]
        if isempty(remaining_tokens) # Should be caught by length check, but defensive
            parse_error("Syntax error near 'not'. Missing operand.")
        end
        # Prevent "not and", "not or", "not not" if "not" is not a general prefix operator in this DSL
        if is_operator(remaining_tokens[begin]) && remaining_tokens[begin] != "not" # allow "not not" if desired, though unusual
             parse_error("Operator '$(remaining_tokens[begin])' cannot directly follow 'not'.")
        end
        return (!, parse_query_vector(remaining_tokens))

    # Base case: No top-level logical operators. Must be a keyword phrase.
    else
        if isempty(s_vec) # Should not happen if initial checks are correct
            parse_error("Unexpected empty query segment.")
        end
        token_keyword_name = s_vec[begin]

        # Standard Keywords (e.g., "name", "resnum", "index")
        for key_obj in keywords # key_obj is of type Keyword
            if token_keyword_name == key_obj.name
                if length(s_vec) == 1 # Keyword name token only, no arguments
                    parse_error("Keyword '$(key_obj.name)' requires at least one argument.")
                end
                
                keyword_args = s_vec[begin+1:end] # Arguments following the keyword name
        
                is_operator_syntax_match = false
                if !isempty(keyword_args)
                    first_arg = keyword_args[1]
                    for op_tuple in key_obj.operators # e.g., ("<", isless)
                        operator_string = op_tuple[1]
                        if first_arg == operator_string
                            # Expected form: "keyword operator value", so keyword_args should be ["operator", "value"] (length 2)
                            if length(keyword_args) == 2
                                is_operator_syntax_match = true
                            else
                                parse_error(
                                    "Malformed operator expression for keyword '$(key_obj.name)'. "*
                                    "Expected 'keyword $operator_string value'. Got: $(join(s_vec, " "))"
                                )
                            end
                            break # Operator string found and processed
                        end
                    end
                end
        
                if is_operator_syntax_match
                    # Case: "keyword operator value", e.g., "resnum < 13"
                    # keyword_args will be ["<", "13"]. The Keyword functor handles this structure.
                    return key_obj(keyword_args)
                else
                    # Case: Not a recognized "keyword operator value" structure.
                    # This implies implicit equality: "keyword value" or "keyword value1 value2 ..." (for OR expansion).
        
                    if isempty(keyword_args) # Should have been caught by length(s_vec) == 1
                         parse_error("No arguments provided for keyword '$(key_obj.name)'.") # Should be unreachable
                    end
        
                    # Sanity check for multi-value: ensure no operators are present in the value list.
                    # E.g. "resnum 10 < 20" is an error here because "10" is not an operator,
                    # but "<" appears later in a context expecting only values.
                    for arg_val in keyword_args
                        for op_tuple in key_obj.operators
                            if arg_val == op_tuple[1] # op_tuple[1] is the operator string
                                parse_error(
                                    "Syntax error for keyword '$(key_obj.name)'. Operator '$(op_tuple[1])' found in an unexpected position. "*
                                    "Arguments: $(join(keyword_args, " ")). Operator expressions must be 'keyword $(op_tuple[1]) value'."
                                )
                            end
                        end
                    end
        
                    # Proceed with implicit equality (single value or multi-value OR).
                    if length(keyword_args) == 1
                        # e.g., "name CA" -> keyword_args = ["CA"]
                        # The Keyword functor handles this as implicit equality.
                        return key_obj(keyword_args)
                    else
                        # Multi-value implicit OR case, e.g., "resname ARG GLU ASP"
                        # keyword_args = ["ARG", "GLU", "ASP"]
                        current_expr_tree = key_obj([keyword_args[end]]) # Process the last value
                        for k_idx in (length(keyword_args)-1):-1:firstindex(keyword_args) # Iterate remaining values
                            current_expr_tree = (|, key_obj([keyword_args[k_idx]]), current_expr_tree)
                        end
                        return current_expr_tree
                    end
                end
            end
        end
        
        # Functional Keywords (e.g., "element")
        for key_obj in functional_keywords # key_obj is of type FunctionalKeyword
            if token_keyword_name == key_obj.name
                if length(s_vec) == 1
                     parse_error("Functional keyword '$(key_obj.name)' requires at least one argument.")
                end
                
                keyword_args = s_vec[begin+1:end]
        
                if isempty(keyword_args) # Should be caught by length(s_vec) == 1
                    parse_error("No arguments provided for functional keyword '$(key_obj.name)'.") # Should be unreachable
                end
                
                # FunctionalKeywords, as defined, don't have an 'operators' field for infix ops.
                # They generally take a list of values.
                # "element C" or "element C N O" (interpreted as OR).
        
                if length(keyword_args) == 1
                    # e.g., key_obj(["C"])
                    return key_obj(keyword_args)
                else
                    # Multi-value OR case, e.g., "element C N O"
                    current_expr_tree = key_obj([keyword_args[end]])
                    for k_idx in (length(keyword_args)-1):-1:firstindex(keyword_args)
                        current_expr_tree = (|, key_obj([keyword_args[k_idx]]), current_expr_tree)
                    end
                    return current_expr_tree
                end
            end
        end

        # Macro Keywords (e.g., "protein", "water")
        for key_obj in macro_keywords
            if token_keyword_name == key_obj.name
                if length(s_vec) > 1
                    parse_error("Macro keyword '$(key_obj.name)' does not take arguments. Unexpected tokens: $(join(s_vec[begin+1:end], " "))")
                end
                # MacroKeyword functor expects an argument list (empty for macros)
                return key_obj(String[]) 
            end
        end
        
        # Functional Keywords (e.g., "element")
        for key_obj in functional_keywords
            if token_keyword_name == key_obj.name
                if length(s_vec) == 1
                     parse_error("Functional keyword '$(key_obj.name)' requires at least one argument.")
                end
                values = s_vec[begin+1:end]

                if length(values) == 1
                    return key_obj(values)
                else # Multiple values, similar to standard Keyword
                    current_expr_tree = key_obj([values[end]])
                    for k in length(values)-1:-1:firstindex(values)
                        current_expr_tree = (|, key_obj([values[k]]), current_expr_tree)
                    end
                    return current_expr_tree
                end
            end
        end
        
        parse_error("Unknown keyword or invalid syntax at: '$(join(s_vec, " "))'")
    end
end

@testitem "Selections" begin

    atoms = read_pdb(PDBTools.TESTPDB)

    @test length(select(atoms, "name CA")) == 104
    sel = select(atoms, "index = 13")
    @test length(sel) == 1
    @test sel[1].index == 13
    @test sel[1].index_pdb == 13

    @test length(select(atoms, "index > 1 and index < 13")) == 11
    @test length(select(atoms, at -> at.index > 1 && at.index < 13)) == 11
    @test length(select(atoms, "protein")) == 1463
    @test length(select(atoms, isprotein)) == 1463
    @test length(select(atoms, "water")) == 58014
    @test length(select(atoms, iswater)) == 58014
    @test length(select(atoms, "resname GLY")) == 84
    @test length(select(atoms, at -> resname(at) == "GLY")) == 84
    @test length(select(atoms, "segname PROT")) == 1463
    @test length(select(atoms, "residue = 2")) == 11
    @test length(select(atoms, "neutral")) == 1233
    @test length(select(atoms, "charged")) == 230
    @test length(select(atoms, "sidechain")) == 854
    @test length(select(atoms, "acidic")) == 162
    @test length(select(atoms, "basic")) == 68
    @test length(select(atoms, "hydrophobic")) == 399
    @test length(select(atoms, "aliphatic")) == 379
    @test length(select(atoms, "aromatic")) == 344
    @test length(select(atoms, "polar")) == 1008
    @test length(select(atoms, "nonpolar")) == 455
    @test maxmin(atoms, "chain A").xlength ≈ [83.083, 83.028, 82.7] atol = 1e-3

    # Advanced selections
    @test length(select(atoms, "name CA and (residue < 15 or residue > 16)")) == 102
    @test length(select(atoms, "(name CA and residue < 15) or (name N and chain A)")) == 299
    @test length(select(atoms, "(not protein and not water) and (resname TMAO or (resname SOD and index 1470))")) == 2535
    @test length(select(atoms, "not protein and not water or (chain A and resnum < 10)")) == 2662
    @test length(select(atoms, "not protein and not water or (chain A and resnum <= 10)")) == 2673
    @test length(select(atoms, "name CA and resname ALA ARG GLU")) == 14
    @test length(select(atoms, "resname ALA ARG GLU and name N")) == 14
    @test length(select(atoms, "(resname ALA ARG GLU) and (name N or name CA)")) == 28
    @test length(select(atoms, "index 2 3 4 5")) == 4
    @test length(select(atoms, "element C N")) == 1331
    @test length(select(atoms, "not protein and element C N")) == 724

    # malformed expression
    @test_throws ArgumentError select(atoms, "name CA and (residue 1")
    @test_throws ArgumentError select(atoms, "index <")
    @test_throws ArgumentError select(atoms, "index < 1.0")
    @test_throws ArgumentError select(atoms, "indes 1")
    @test_throws ArgumentError select(atoms, "element")
    @test_throws ArgumentError select(atoms, "index 1 element")
    @test_throws ArgumentError select(atoms, "protein 1")
    @test_throws ArgumentError select(atoms, "protein = 1")
    @test_throws ArgumentError select(atoms, "residue 1 < 5")
    @test_throws ArgumentError select(atoms, "residue A")
    @test_throws ArgumentError select(atoms, "residue 1 and ()")

    # test string interpolation
    t = "protein"
    n = "CA"
    r = 15
    at = filter(sel"$t and name $n and residue $r", atoms)
    @test length(at) == 1
    @test name(at[1]) == "CA"
    @test residue(at[1]) == 15
    @test_throws ArgumentError filter(sel"$t and name $n and residue abc", atoms)
    @test_throws UndefVarError filter(sel"$t and name $n and residue $r and $MM", atoms)

    # Test editing a field
    atoms[1].index = 0
    @test atoms[1].index == 0

    # Test residue iterator
    someresidues = select(atoms, "residue < 15")
    let
        n = 0
        m = 0.0
        for res in eachresidue(someresidues)
            if name(res) == "SER"
                n += 1
                for atom in res
                    m += mass(atom)
                end
            end
        end
        @test n == 4
        @test m ≈ 348.31340000000006
    end

    # Residue properties (discontinous set)
    lessresidues = select(someresidues, "residue < 3 or residue > 12")
    @test residue.(eachresidue(lessresidues)) == [1, 2, 13, 14]
    @test resnum.(eachresidue(lessresidues)) == [1, 2, 13, 14]
    @test name.(eachresidue(lessresidues)) == ["ALA", "CYS", "SER", "SER"]
    @test resname.(eachresidue(lessresidues)) == ["ALA", "CYS", "SER", "SER"]
    @test chain.(eachresidue(lessresidues)) == ["A", "A", "A", "A"]
    @test model.(eachresidue(lessresidues)) == [1, 1, 1, 1]

end # testitem