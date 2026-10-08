"""
    VolumetricData(data::AbstractArray{<:Real,3}; origin, step)

Values of a property on a regular three-dimensional grid (for example, a density map).
`data[i, j, k]` is the value at position `origin .+ step .* (i - 1, j - 1, k - 1)`,
where `origin` is the position of the first grid point, and `step` the spacing of the
grid in each direction (a number, or a vector of three numbers).

Volumetric data can be read from and written to OpenDX (`.dx`) files with
[`read_dx`](@ref) and [`write_dx`](@ref), and displayed as isosurfaces with [`visualize`](@ref).

# Example

```jldoctest
julia> using PDBTools

julia> v = VolumetricData(zeros(10, 20, 30); origin=[0.0, 0.0, 0.0], step=0.5)
PDBTools.VolumetricData: 10×20×30 grid, origin = [0.0, 0.0, 0.0], step = [0.5, 0.5, 0.5]
    values from 0.0 to 0.0
```

"""
struct VolumetricData
    data::Array{Float64,3}
    origin::SVector{3,Float64}
    step::SVector{3,Float64}
end
function VolumetricData(data::AbstractArray{<:Real,3}; origin, step)
    length(origin) == 3 || throw(ArgumentError("origin must have three components. Got: $origin"))
    step = step isa Real ? SVector(step, step, step) : step
    length(step) == 3 || throw(ArgumentError("step must be a number or have three components. Got: $step"))
    Base.all(>(0), step) || throw(ArgumentError("step must be positive. Got: $step"))
    return VolumetricData(Array{Float64,3}(data), SVector{3,Float64}(origin), SVector{3,Float64}(step))
end

function Base.show(io::IO, ::MIME"text/plain", v::VolumetricData)
    print(io, chomp("""
    PDBTools.VolumetricData: $(join(size(v.data), "×")) grid, origin = $(Vector(v.origin)), step = $(Vector(v.step))
        values from $(minimum(v.data)) to $(maximum(v.data))
    """))
end

"""
    read_dx(filename::AbstractString)

Reads volumetric data from a file in the OpenDX (`.dx`) format, as written by APBS, VMD,
and other programs. Returns a [`VolumetricData`](@ref) object. Only orthogonal grids are supported.

"""
function read_dx(filename::AbstractString)
    n = nothing
    origin = nothing
    deltas = SVector{3,Float64}[]
    values = Float64[]
    open(expanduser(filename)) do io
        reading_data = false
        for line in eachline(io)
            line = strip(line)
            (isempty(line) || startswith(line, "#")) && continue
            if reading_data
                if startswith(line, "attribute") || startswith(line, "object") || startswith(line, "component")
                    reading_data = false
                    continue
                end
                for v in split(line)
                    push!(values, parse(Float64, v))
                end
            elseif startswith(line, "object") && occursin("gridpositions", line)
                n = Tuple(parse.(Int, split(line)[end-2:end]))
            elseif startswith(line, "origin")
                origin = SVector{3,Float64}(parse.(Float64, split(line)[2:4]))
            elseif startswith(line, "delta")
                push!(deltas, SVector{3,Float64}(parse.(Float64, split(line)[2:4])))
            elseif occursin("data follows", line)
                reading_data = true
            end
        end
    end
    if isnothing(n) || isnothing(origin) || length(deltas) != 3
        throw(ArgumentError("Could not read the grid definition from $filename."))
    end
    for i in 1:3, j in 1:3
        i != j && deltas[i][j] != 0 && throw(ArgumentError("Only orthogonal grids are supported."))
    end
    length(values) == prod(n) ||
        throw(ArgumentError("Number of values ($(length(values))) does not match the grid size $(n) in $filename."))
    # In the DX format, the last index varies fastest
    data = permutedims(reshape(values, reverse(n)), (3, 2, 1))
    return VolumetricData(data, origin, SVector(deltas[1][1], deltas[2][2], deltas[3][3]))
end

"""
    write_dx(filename::AbstractString, v::VolumetricData)

Writes the volumetric data `v` to a file in the OpenDX (`.dx`) format, which can be read
by visualization software (VMD, PyMOL, ChimeraX, 3Dmol.js).

"""
function write_dx(filename::AbstractString, v::VolumetricData)
    open(expanduser(filename), "w") do io
        _write_dx(io, v)
    end
    return filename
end

function _write_dx(io::IO, v::VolumetricData)
    n = size(v.data)
    println(io, "# Created by PDBTools.jl")
    println(io, "object 1 class gridpositions counts $(n[1]) $(n[2]) $(n[3])")
    @printf(io, "origin %.6f %.6f %.6f\n", v.origin...)
    @printf(io, "delta %.6f 0 0\n", v.step[1])
    @printf(io, "delta 0 %.6f 0\n", v.step[2])
    @printf(io, "delta 0 0 %.6f\n", v.step[3])
    println(io, "object 2 class gridconnections counts $(n[1]) $(n[2]) $(n[3])")
    println(io, "object 3 class array type double rank 0 items $(prod(n)) data follows")
    # The last index varies fastest
    icol = 0
    for i in 1:n[1], j in 1:n[2], k in 1:n[3]
        x = v.data[i, j, k]
        print(io, x == 0 ? "0" : @sprintf("%.6g", x))
        icol += 1
        print(io, icol % 3 == 0 ? "\n" : " ")
    end
    icol % 3 == 0 || println(io)
    println(io, "attribute \"dep\" string \"positions\"")
    println(io, "object \"density\" class field")
    println(io, "component \"positions\" value 1")
    println(io, "component \"connections\" value 2")
    println(io, "component \"data\" value 3")
    return io
end

@testitem "VolumetricData and DX files" begin
    using PDBTools
    data = reshape(Float64.(1:24), 2, 3, 4)
    data[1, 1, 1] = 0.0
    v = VolumetricData(data; origin=[1.0, 2.0, 3.0], step=0.5)
    @test v.step == [0.5, 0.5, 0.5]
    @test occursin("2×3×4 grid", sprint(show, MIME"text/plain"(), v))
    file = tempname() * ".dx"
    @test write_dx(file, v) == file
    lines = readlines(file)
    @test lines[2] == "object 1 class gridpositions counts 2 3 4"
    # last index varies fastest
    @test split(lines[9]) == ["0", "7", "13"]
    v2 = read_dx(file)
    @test v2.data == v.data
    @test v2.origin == v.origin
    @test v2.step == v.step
    v = VolumetricData(data; origin=(0, 0, 0), step=[0.5, 1.0, 2.0])
    write_dx(file, v)
    @test read_dx(file).step == [0.5, 1.0, 2.0]
    @test_throws ArgumentError VolumetricData(data; origin=[0.0, 0.0], step=1.0)
    @test_throws ArgumentError VolumetricData(data; origin=[0.0, 0.0, 0.0], step=-1.0)
    @test_throws ArgumentError VolumetricData(data; origin=[0.0, 0.0, 0.0], step=[1.0, 1.0])
    write(file, "object 1 class gridpositions counts 2 2 2\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 3 data follows\n1 2 3\n")
    @test_throws ArgumentError read_dx(file)
    write(file, "object 1 class gridpositions counts 1 1 1\norigin 0 0 0\ndelta 1 1 0\ndelta 0 1 0\ndelta 0 0 1\nobject 3 data follows\n1\n")
    @test_throws ArgumentError read_dx(file)
    write(file, "origin 0 0 0\n")
    @test_throws ArgumentError read_dx(file)
    rm(file)
end
