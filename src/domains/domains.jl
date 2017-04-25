# domains.jl

"An N-dimensional domain."
abstract AbstractDomain{N}

ndims{N}(::Type{AbstractDomain{N}}) = N
ndims{D <: AbstractDomain}(::Type{D}) = ndims(supertype(D))
ndims{N}(::AbstractDomain{N}) = N


typealias AbstractDomain1d AbstractDomain{1}
typealias AbstractDomain2d AbstractDomain{2}
typealias AbstractDomain3d AbstractDomain{3}
typealias AbstractDomain4d AbstractDomain{4}

# left and right of domains falls back to bounding box domains
left(d::AbstractDomain) = left(boundingbox(d))
right(d::AbstractDomain) = right(boundingbox(d))

left(d::AbstractDomain, i::Int) = left(boundingbox(d), i)
right(d::AbstractDomain, i::Int) = right(boundingbox(d), i)

# Redirect calls to 'in' to 'indomain', because the latter has fewer
# methods than the standard function 'in' and that makes it easier to do
# duck typing.
in{N}(x::SVector{N}, d::AbstractDomain{N}) = indomain(x, d)
in(x::Number, d::AbstractDomain1d) = indomain(x, d)
in(x::SVector{1}, d::AbstractDomain1d) = indomain(x[1], d)

# Convert a point given as any other vector into a SVector
in{N}(x::AbstractVector, d::AbstractDomain{N}) = in(SVector{N}(x), d)

# Check whether a value is in an interval, up to 10 times machine precision
in{T <: AbstractFloat}(x::Number, a::T, b::T) = (a-10eps(T) <= x <= b+10eps(T))
in{T <: Number}(x::Number, a::T, b::T) = a <= x <= b

# Evaluation on a grid should be implemented by indomain_grid for each domain
in{N}(g::AbstractGrid{N}, d::AbstractDomain{N}) = indomain_grid(g, d)

# Default methods for evaluation on a grid: the default is to call eval on the domain with
# points as arguments. Domains that have faster grid evaluation routines may define their own version.
indomain_grid(g::AbstractGrid, d::AbstractDomain) = indomain_grid!(zeros(Bool, size(g)), g, d)

# Note that indomain_grid! only updates the result - it should be initialized to all false!
# The idea is that you can chain different calls to indomain_grid (as used in DomainCollection below)
function indomain_grid!{N}(result, grid::AbstractGrid{N}, domain::AbstractDomain{N})
    for (i,x) in enumerate(grid)
        result[i] |= indomain(x, domain)
    end
    result
end


## Arithmetics

# Make sure domains only need to implement addition/multiplication with numbers to the right
(+)(x::Number, d::AbstractDomain) = d + x
(+)(x::AbstractVector, d::AbstractDomain) = d + x
(*)(x::Number, d::AbstractDomain) = d * x

(/)(d::AbstractDomain, x::Number) = d * (1/x)


include("tensorproductdomain.jl")
include("derived_domains.jl")
include("specific_domains.jl")
include("curve.jl")
