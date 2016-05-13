# domainframe.jl

"""
A DomainFrame is the restriction of a basis to a subset of its domain. This results
in a frame.
"""
immutable DomainFrame{N,T} <: AbstractFrame{N,T}
    domain      ::  AbstractDomain{N}
    basis       ::  FunctionSet{N,T}

    function DomainFrame(domain::AbstractDomain, basis::FunctionSet)
        @assert is_basis(basis)

        new(domain, basis)
    end
end

DomainFrame{N,T}(domain::AbstractDomain{N}, basis::FunctionSet{N,T}) =
    DomainFrame{N,T}(domain, basis)

basis(f::DomainFrame) = f.basis

domain(f::DomainFrame) = f.domain

name(f::DomainFrame) = "A frame of " * name(f.basis)

promote_eltype{S}(f::DomainFrame, ::Type{S}) =
    DomainFrame(f.domain, promote_eltype(f.basis, S))

resize(f::DomainFrame, n) = DomainFrame(domain(f), resize(basis(f), n))

extension_size(f::DomainFrame) = extension_size(basis(f))


for op in (:size, :length)
    @eval $op(f::DomainFrame, args...) = $op(f.basis, args...)
end

for op in (:differentiation_operator, :antidifferentiation_operator)
    @eval function $op(f::DomainFrame, order)
        op = $op(basis(f), order)
        WrappedOperator(f,DomainFrame(domain(f),dest(op)),op)
    end
end

for op in (:extension_operator,)
    @eval function $op(f1::DomainFrame, f2::DomainFrame)
        @assert is_compatible(f1,f2)
        op = $op(basis(f1),basis(f2))
        WrappedOperator(f1,f2,op)
    end
end


is_compatible(d1::DomainFrame, d2::DomainFrame) = is_compatible(basis(d1),basis(d2))

function (*)(d1::DomainFrame, d2::DomainFrame, args...)
    @assert is_compatible(d1,d2) 
    (mset, mcoef) = (*)(basis(d1),basis(d2),args...)
    df = DomainFrame(domain(d1) ∩ domain(d2),mset)
    (df, mcoef)
end

# Should we check whether x lies in the domain?
call_set(e::SetExpansion, s::DomainFrame, coef, x...) = call_expansion(basis(s), coef, x...)

call_set!(result, e::SetExpansion, s::DomainFrame, coef, x...) = call_expansion!(result, basis(s), coef, x...)

call_element(s::DomainFrame, idx::Int, x...) = call_element(basis(s), idx, x...)

"""
Make a DomainFrame, but match tensor product domains with tensor product sets in a suitable way.

For example: an interval ⊗ a disk (= a cylinder) combined with a 3D Fourier series, leads to a
tensor product of a Fourier series on the interval ⊗ a 2D Fourier series on the disk.
"""
function domainframe(domain::TensorProductDomain, basis::TensorProductSet)
    domainframes = FunctionSet[]
    dc = 1
    for i = 1:composite_length(domain)
        el = element(domain, i)
        range = dc:dc+dim(el)-1
        push!(domainframes, DomainFrame(el, element(basis, range)))
        dc += dim(el)
    end
    tensorproduct(domainframes...)
end

domainframe(domain::AbstractDomain, basis::FunctionSet) = DomainFrame(domain, basis)
