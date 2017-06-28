# tensorproductdomain.jl

###############################################################################################
### A Tensor Product of Domains
###############################################################################################

"""
A TensorProductDomain represents the tensor product of other domains.

immutable TensorProductDomain{TD,N} <: AbstractDomain{N}

Parameters:
- TD is a tuple of (domain) types
- N is the total dimension of this domain
"""
immutable TensorProductDomain{TD,N} <: AbstractDomain{N}
	domains	::	TD
end

# Generic functions for composite types:
elements(d::TensorProductDomain) = d.domains
element(d::TensorProductDomain, j::Int) = d.domains[j]
composite_length(d::TensorProductDomain) = length(d.domains)

function TensorProductDomain(domains...)
    TD = typeof(domains)
    N = sum(map(ndims, domains))
    TensorProductDomain{TD,N}(domains)
end

tensorproduct(d::AbstractDomain) = d
tensorproduct(d::AbstractDomain, n::Int) = tensorproduct([d for i=1:n]...)
tensorproduct(d::AbstractDomain...) =
    TensorProductDomain(flatten(TensorProductDomain, d...)...)

^(d::AbstractDomain, n::Int) = tensorproduct(d, n)

indomain(x, t::TensorProductDomain) = indomain(x, elements(t)...)

indomain(x::SVector{2}, d1::AbstractDomain{1}, d2::AbstractDomain{1}) =
	indomain(x[1], d1) && indomain(x[2], d2)

indomain(x::SVector{3}, d1::AbstractDomain{1}, d2::AbstractDomain{1}, d3::AbstractDomain{1}) =
	indomain(x[1], d1) && indomain(x[2], d2) && indomain(x[3], d3)

indomain(x::SVector{4}, d1::AbstractDomain{1}, d2::AbstractDomain{1}, d3::AbstractDomain{1}, d4::AbstractDomain{1}) =
	indomain(x[1], d1) && indomain(x[2], d2) && indomain(x[3], d3) && indomain(x[4], d4)

indomain(x::SVector{3}, d1::AbstractDomain{1}, d2::AbstractDomain{2}) =
	indomain(x[1], d1) && indomain(SVector(x[2],x[3]), d2)

indomain(x::SVector{3}, d1::AbstractDomain{2}, d2::AbstractDomain{1}) =
	indomain(SVector(x[1],x[2]), d1) && indomain(x[3], d2)

indomain(x::SVector{4}, d1::AbstractDomain{2}, d2::AbstractDomain{2}) =
	indomain(SVector(x[1],x[2]), d1) && indomain(SVector(x[3],x[4]), d2)

indomain(x::SVector{4}, d1::AbstractDomain{1}, d2::AbstractDomain{3}) =
	indomain(x[1], d1) && indomain(SVector(x[2],x[3],x[4]), d2)

indomain(x::SVector{4}, d1::AbstractDomain{3}, d2::AbstractDomain{1}) =
	indomain(SVector(x[1],x[2],x[3]), d1) && indomain(x[1], d2)

indomain(x::SVector{4}, d1::AbstractDomain{1}, d2::AbstractDomain{1}, d3::AbstractDomain{2}) =
	indomain(x[1], d1) && indomain(x[2], d2) && indomain(SVector(x[3],x[4]), d3)

indomain(x::SVector{4}, d1::AbstractDomain{1}, d2::AbstractDomain{2}, d3::AbstractDomain{1}) =
	indomain(x[1], d1) && indomain(SVector(x[2],x[3]), d2) && indomain(x[4], d3)

indomain(x::SVector{4}, d1::AbstractDomain{2}, d2::AbstractDomain{1}, d3::AbstractDomain{1}) =
	indomain(SVector(x[1],x[2]), d1) && indomain(x[3], d2) && indomain(x[4], d3)

# TODO: make this code for indomain more general!

# TODO: provide implementation of indomain for tensorproductgrids

boundingbox(d::TensorProductDomain) = tensorproduct(map(boundingbox, elements(d))...)

function show(io::IO, t::TensorProductDomain)
    L = composite_length(t)
    for i in 1:L-1
        show(io, element(t, i))
        print(io, " x ")
    end
    show(io, element(t, L))
end


(*)(d::TensorProductDomain, x::Number) = tensorproduct([domain*x for domain in elements(d)]...)

dist(x, t::TensorProductDomain) = minimum(map(dist,x,elements(t)))

function normal(x, t::TensorProductDomain)
    index = findmin(map(dist,x,elements(t)))[2]
    [(i==index)*normal(x[i],element(t,i)) for i =1:composite_length(t)]
end
