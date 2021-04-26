#
# Function to return a three-letter code residue name from the one-letter or residue names
#
function threeletter(resname::String)
  if length(resname) == 1
    ires = findfirst(r->lowercase(r.one_letter_code) == lowercase(resname), natural_aminoacids)
  else
    ires = findfirst(r->lowercase(r.name) == lowercase(resname), natural_aminoacids)
  end
  if ires == nothing
    return nothing
  else
    return natural_aminoacids[ires].three_letter_code
  end
end

