using PDBTools
using Test

tests = [ "./select.jl" ]

for test in tests
  include(test)
end

