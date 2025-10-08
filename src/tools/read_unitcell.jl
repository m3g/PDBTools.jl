"""
    read_unitcell(file::AbstractString)

Reads the lattice parameters of unitcell from a PDB (`CRYST1` field) or mmCIF file, and converts it to a unitcell matrix
with the `lattice_to_matrix` function.

# Example

```julia-repl
julia> using PDBTools

julia> m = read_unitcell(PDBTools.TESTPBC)
3×3 StaticArraysCore.SMatrix{3, 3, Float32, 9} with indices SOneTo(3)×SOneTo(3):
 85.0  -3.71547f-6  -3.71547f-6
  0.0  85.0         -3.71547f-6
  0.0   0.0         85.

julia> matrix_to_lattice(m)
(a = 85.0f0, b = 85.0f0, c = 85.0f0, α = 90.0f0, β = 90.0f0, γ = 90.0f0)
```

"""
@views function read_unitcell(file::AbstractString)
    a = nothing
    b = nothing 
    c = nothing
    α = nothing
    β = nothing
    γ = nothing
    open(file) do io
        for line in eachline(io)
            ifirst = findfirst(!isspace, line)
            ilast = findlast(!isspace, line)
            if ilast - ifirst + 1 >= 6
                key = line[ifirst:ifirst+5]
                # File is a PDB file
                if key == "CRYST1"
                    uc = split(line)
                    a = parse(Float32, uc[2])
                    b = parse(Float32, uc[3])
                    c = parse(Float32, uc[4])
                    α = parse(Float32, uc[5])
                    β = parse(Float32, uc[6])
                    γ = parse(Float32, uc[7])
                    break
                end
                # mmCIF file
                if key == "_cell."
                    inext_space = findnext(isspace, line, ifirst)
                    val = line[ifirst+6:inext_space-1]
                    val == "length_a" && (a = parse(Float32, line[inext_space+1:end]))
                    val == "length_b" && (b = parse(Float32, line[inext_space+1:end]))
                    val == "length_c" && (c = parse(Float32, line[inext_space+1:end]))
                    val == "angle_alpha" && (α = parse(Float32, line[inext_space+1:end]))
                    val == "angle_beta" && (β = parse(Float32, line[inext_space+1:end]))
                    val == "angle_gamma" && (γ = parse(Float32, line[inext_space+1:end]))
                end
                Base.all(!isnothing, (a, b, c, α, β, γ)) && break
            end
        end
    end
    if any(isnothing, (a, b, c, α, β, γ))
        throw(ArgumentError("Could not find (complete?) unit cell information in file."))
    end
    return lattice_to_matrix(a, b, c, α, β, γ)
end


"""
    lattice_to_matrix(a, b, c, α, β, γ)

Converts unit cell lattice parameters and convert them to a 3x3
unit cell matrix (orthogonalization matrix).

The resulting matrix has the lattice vectors as its columns, with
vector 'a' aligned along the x-axis and vector 'b' in the xy-plane.
This matrix can be used to transform fractional coordinates to
Cartesian coordinates.

# Arguments
- `a::Real`: Length of side a in Ångströms.
- `b::Real`: Length of side b in Ångströms.
- `c::Real`: Length of side c in Ångströms.
- `α::Real`: Angle alpha in degrees.
- `β::Real`: Angle beta in degrees.
- `γ::Real`: Angle gamma in degrees.

# Returns
- A 3x3 static matrix representing the unit cell vectors.

"""
function lattice_to_matrix(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    a, b, c, α, β, γ = promote(a, b, c, α, β, γ)

    # Convert angles from degrees to radians for trigonometric functions
    α_rad, β_rad, γ_rad = deg2rad(α), deg2rad(β), deg2rad(γ)

    cos_α, cos_β, cos_γ = cos(α_rad), cos(β_rad), cos(γ_rad)
    sin_γ = sin(γ_rad)

    # Calculate the unit cell volume V = abc * ω
    # ω is a term related to the Gram determinant
    ω = sqrt(1 - cos_α^2 - cos_β^2 - cos_γ^2 + 2 * cos_α * cos_β * cos_γ)

    # Construct the transformation matrix M
    T = typeof(a)
    M = SMatrix{3,3,T,9}(
        a, 
        zero(T),
        zero(T),
        b * cos_γ,
        b * sin_γ,
        zero(T),
        c * cos_β,
        c * (cos_α - cos_β * cos_γ) / sin_γ,
        c * ω / sin_γ
    )

    return M
end

"""
    matrix_to_lattice(M::AbstractMatrix)

Converts a 3x3 unit cell matrix back into the six lattice parameters 
(sides a, b, c and angles α, β, γ).

# Arguments
- `M::AbstractMatrix`: A 3x3 matrix where columns are the lattice vectors.

# Returns
- `NamedTuple`: A named tuple containing the six parameters:
  `(a, b, c, α, β, γ)`. Lengths are in the matrix's original units
  (typically Ångströms), and angles are in degrees.

"""
function matrix_to_lattice(M::AbstractMatrix)
    T = eltype(M)
    # Ensure the input is a 3x3 matrix
    if size(M) != (3, 3)
        throw(ArgumentError("Input must be a 3x3 matrix."))
    end

    # Extract column vectors
    a_vec = @view(M[:, 1])
    b_vec = @view(M[:, 2])
    c_vec = @view(M[:, 3])

    # Calculate side lengths
    a = norm(a_vec)
    b = norm(b_vec)
    c = norm(c_vec)

    # Calculate angles in degrees
    # We use clamp() to avoid domain errors with acos due to floating point inaccuracies
    cos_gamma = clamp(dot(a_vec, b_vec) / (a * b), -one(T), one(T))
    cos_beta  = clamp(dot(a_vec, c_vec) / (a * c), -one(T), one(T))
    cos_alpha = clamp(dot(b_vec, c_vec) / (b * c), -one(T), one(T))

    α = rad2deg(acos(cos_alpha))
    β = rad2deg(acos(cos_beta))
    γ = rad2deg(acos(cos_gamma))

    return (a=a, b=b, c=c, α=α, β=β, γ=γ)
end

@testitem "read_unitcell" begin
    using PDBTools

    a, b, c = 54.530, 63.780, 82.340
    α, β, γ = 90.0, 90.0, 90.0
    ortho_matrix = lattice_to_matrix(a, b, c, α, β, γ)
    @test ortho_matrix ≈ [54.530 0 0; 0 63.780 0; 0 0 82.340]
    @test eltype(ortho_matrix) == Float64
    @test all(v1 ≈ v2 for (v1,v2) in zip(values(matrix_to_lattice(ortho_matrix)),values((a = 54.53, b = 63.78, c = 82.34, α = 90.0, β = 90.0, γ = 90.0))))
    a, b, c = 54.530f0, 63.780f0, 82.340f0
    α, β, γ = 90.0f0, 90.0f0, 90.0f0
    ortho_matrix = lattice_to_matrix(a, b, c, α, β, γ)
    @test eltype(ortho_matrix) == Float32

    a, b, c = 39.850, 54.450, 59.840
    α, β, γ = 72.01, 89.44, 89.43
    triclinic_matrix = lattice_to_matrix(a, b, c, α, β, γ)
    @test triclinic_matrix ≈ [39.85 0.54168 0.584858; 0 54.4473 18.4767; 0 0 56.913] atol=1e-4
    @test eltype(triclinic_matrix) == Float64
    @test all(isapprox(v1,v2,atol=1e-5) for (v1,v2) in zip(values(matrix_to_lattice(triclinic_matrix)),values((a = 39.85, b = 54.449999999999996, c = 59.84, α = 72.01, β = 89.44, γ = 89.43))))

    @test read_unitcell(PDBTools.TESTPBC) ≈ [107.845 0 0;  0.0  107.845 0; 0.0 0.0 107.845]
    @test read_unitcell(PDBTools.TESTCIF) ≈ [74.449 -3.25427f-6 -3.25427f-6; 0.0 74.449 -3.25427f-6; 0.0 0.0 74.449]
    
    # Throw error if unitcel ldata is missing
    @test_throws ArgumentError read_unitcell(PDBTools.TESTPDB)
    @test_throws ArgumentError read_unitcell(PDBTools.BROKENCIF)
    @test_throws ArgumentError matrix_to_lattice([1 0; 0 1])

    # Test allocations of conversions
    using Chairmarks
    a, b, c = 39.850, 54.450, 59.840
    α, β, γ = 72.01, 89.44, 89.43
    t = @b lattice_to_matrix($a, $b, $c, $α, $β, $γ)
    @test t.allocs == 0
    l = lattice_to_matrix(a, b, c, α, β, γ)
    t = @b matrix_to_lattice($l)
    @test t.allocs == 0

end