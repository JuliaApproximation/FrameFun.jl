# domains.jl

"An N-dimensional domain."
abstract AbstractDomain{N}

dim{N}(d::AbstractDomain{N}) = N


typealias AbstractDomain1d AbstractDomain{1}
typealias AbstractDomain2d AbstractDomain{2}
typealias AbstractDomain3d AbstractDomain{3}
typealias AbstractDomain4d AbstractDomain{4}

# left and right of domains falls back to bounding box domains
left(d::AbstractDomain) = left(boundingbox(d))
right(d::AbstractDomain) = right(boundingbox(d))

left(d::AbstractDomain, i::Int) = left(boundingbox(d), i)
right(d::AbstractDomain, i::Int) = right(boundingbox(d), i)

# Domains are evaluated using vectors to specify the points, except in 1D
# Provide fallback routine for users not using vectors in 1d
in(x::Number, d::AbstractDomain1d) = in(Vec(x), d)

# Convert a point given as an array into a Vec point
in{N}(x::AbstractVector, d::AbstractDomain{N}) = in(Vec{N,eltype(x)}(x...), d)

# Check whether a value is in an interval, up to 10 times machine precision
in{T <: AbstractFloat}(x::Number, a::T, b::T) = (a-10eps(T) <= x <= b+10eps(T))
in{T <: Number}(x::Number, a::T, b::T) = a <= x <= b

# Evaluation on a grid should be implemented by evalgrid for each domain
in{N}(g::AbstractGrid{N}, d::AbstractDomain{N}) = evalgrid(g, d)

# Default methods for evaluation on a grid: the default is to call eval on the domain with
# points as arguments. Domains that have faster grid evaluation routines may define their own version.
evalgrid(g::AbstractGrid, d::AbstractDomain) = evalgrid!(zeros(Bool, size(g)), g, d)

# Note that evalgrid! only updates the result - it should be initialized to all false!
# The idea is that you can chain different calls to evalgrid (as used in DomainCollection below)
function evalgrid!{N}(result, g::AbstractGrid{N}, d::AbstractDomain{N})
    for i in eachindex(g)
        result[i] |= in(g[i], d)
    end
    result
end

show{N}(io::IO, v::Vec{N}) = print(io, Vector(v))


## Arithmetics

# Make suredomains only need to implement addition/multiplication with numbers to the right
(+)(x::Number, d::AbstractDomain) = d + x
(+)(x::AnyVector, d::AbstractDomain) = d + x

(*)(x::Number, d::AbstractDomain) = d * x

(/)(d::AbstractDomain, x::Number) = d * (1/x)


include("tensorproductdomain.jl")
include("derived_domains.jl")
include("specific_domains.jl")
