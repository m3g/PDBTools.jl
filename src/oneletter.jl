#
# Function to return a oneletter residue name from the three letter name or residue name
#
function oneletter(resname::String)
  if length(resname) == 3
    ires = findfirst(r->lowercase(r.three_letter_code) == lowercase(resname), natural_aminoacids)
  else
    ires = findfirst(r->lowercase(r.name) == lowercase(resname), natural_aminoacids)
  end
  if ires == nothing
    return nothing
  else
    return natural_aminoacids[ires].one_letter_code
  end
end

