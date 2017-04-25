# randomgrid.jl

"""
Compute a scattered grid of M points randomly distributed in `Ω`, using
the uniform probability measure on `Ω`.
"""
randomgrid{N}(Ω::AbstractDomain{N}, M::Int, T = Float64) =
    ScatteredGrid([randompoint(Ω, boundingbox(Ω), T) for m in 1:M])

"Generate a single random point inside the given box, with `eltype` `T`."
function randompoint{N}(box::BBox{N}, T = Float64)
    vals = Point{N,T}(rand(N))
    vals .* left(box) + (1-vals) .* right(box)
end

# Don't return an SVector in 1d, just a value
function randompoint(box::BBox{1}, T = Float64)
    val = rand()
    val * left(box) + (1-val) * right(box)
end

"""
Generate a single random point inside the given domain, with `eltype` `T`.
Random points are generated inside the given box, until one is inside the domain.
"""
function randompoint{N}(Ω::AbstractDomain{N}, box = boundingbox(Ω), T = Float64)
    local point
    indomain = false
    while ~indomain
        point = randompoint(box, T)
        indomain = point ∈ Ω
    end
    point
end
