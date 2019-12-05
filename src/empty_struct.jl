#
# Returns an empty structure
#

function empty_struct(X :: DataType)
  list = []
  for type in X.types
    if type == Int64
      push!(list,0)
    end
    if type == String
      push!(list,"X")
    end
    if type == Float64
      push!(list,0.)
    end
  end
  y = X(list...)
  return y
end
