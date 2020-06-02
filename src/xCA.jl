#
# Return the coordinates of the CA atoms only in a vector
#

xCA( atoms :: Vector{Atom} ) = xName( atoms, "CA" )

