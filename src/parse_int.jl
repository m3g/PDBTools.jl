#
# Function that tries to read a number as an integer given a number
# or, perhaps, an hexadecimal representation string
#

function parse_int(s :: String)
  try
    i = parse(Int64,s)
    return i
  catch
    i = parse(Int64,s,base=16)
    return i
  end
  error("Could not read integer from string: \"$s\"")
end

