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


in{TD,N}(x::SVector{N}, t::TensorProductDomain{TD,N}) = in(x, elements(t)...)

in(x::SVector{2}, d1::AbstractDomain{1}, d2::AbstractDomain{1}) =
	in(x[1], d1) && in(x[2], d2)

in(x::SVector{3}, d1::AbstractDomain{1}, d2::AbstractDomain{1}, d3::AbstractDomain{1}) =
	in(x[1], d1) && in(x[2], d2) && in(x[3], d3)

in(x::SVector{4}, d1::AbstractDomain{1}, d2::AbstractDomain{1}, d3::AbstractDomain{1}, d4::AbstractDomain{1}) =
	in(x[1], d1) && in(x[2], d2) && in(x[3], d3) && in(x[4], d4)

in(x::SVector{3}, d1::AbstractDomain{1}, d2::AbstractDomain{2}) =
	in(x[1], d1) && in(SVector(x[2],x[3]), d2)

in(x::SVector{3}, d1::AbstractDomain{2}, d2::AbstractDomain{1}) =
	in(SVector(x[1],x[2]), d1) && in(x[3], d2)

in(x::SVector{4}, d1::AbstractDomain{2}, d2::AbstractDomain{2}) =
	in(SVector(x[1],x[2]), d1) && in(SVector(x[3],x[4]), d2)

in(x::SVector{4}, d1::AbstractDomain{1}, d2::AbstractDomain{3}) =
	in(x[1], d1) && in(SVector(x[2],x[3],x[4]), d2)

in(x::SVector{4}, d1::AbstractDomain{3}, d2::AbstractDomain{1}) =
	in(SVector(x[1],x[2],x[3]), d1) && in(x[1], d2)

in(x::SVector{4}, d1::AbstractDomain{1}, d2::AbstractDomain{1}, d3::AbstractDomain{2}) =
	in(x[1], d1) && in(x[2], d2) && in(SVector(x[3],x[4]), d3)

in(x::SVector{4}, d1::AbstractDomain{1}, d2::AbstractDomain{2}, d3::AbstractDomain{1}) =
	in(x[1], d1) && in(SVector(x[2],x[3]), d2) && in(x[4], d3)

in(x::SVector{4}, d1::AbstractDomain{2}, d2::AbstractDomain{1}, d3::AbstractDomain{1}) =
	in(SVector(x[1],x[2]), d1) && in(x[3], d2) && in(x[4], d3)

# # TODO: make this code for in more general!
# # The problem is you can't slice a Vec, so the implementation below for AbstractArray does not work
# # for FixedSizeArray's.
# # All implementations below allocate memory (arrays) except the first one.
# function in{TD,DN,N}(x::SVector{N}, t::TensorProductDomain{TD,DN,N,N})
#     z1 = true
#     for i = 1:N
#         z1 = z1 & in(x[i], t.domains[i])
#     end
#     z1
# end
#
# function in{TD,DN,T}(x::SVector{3,T}, t::TensorProductDomain{TD,DN,2,3})
#     N1 = DN[1]
#     N2 = DN[2]
#     d1 = element(t, 1)
#     d2 = element(t, 2)
#     x1 = SVector{N1,T}([x[j] for j=1:N1])
#     x2 = SVector{N2,T}([x[j] for j=N1+1:N1+N2])
#     in(x1, d1) && in(x2, d2)
# end
#
# function in{TD,DN,T}(x::SVector{4,T}, t::TensorProductDomain{TD,DN,2,4})
#     N1 = DN[1]
#     N2 = DN[2]
#     d1 = element(t, 1)
#     d2 = element(t, 2)
#     x1 = SVector{N1,T}([x[j] for j=1:N1])
#     x2 = SVector{N2,T}([x[j] for j=N1+1:N1+N2])
#     in(x1, d1) && in(x2, d2)
# end
#
# function in{TD,DN,T}(x::SVector{4,T}, t::TensorProductDomain{TD,DN,3,4})
#     d1 = element(t, 1)
#     d2 = element(t, 2)
#     d3 = element(t, 3)
#     x1 = SVector{N1,T}([x[j] for j=1:N1])
#     x2 = SVector{N2,T}([x[j] for j=N1+1:N1+N2])
#     x3 = SVector{N3,T}([x[j] for j=N1+N2+1:N1+N2+N3])
#     in(x1, d1) && in(x2, d2) && in(x3, d3)
# end


# TODO: provide implementation of in for tensorproductgrids

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
