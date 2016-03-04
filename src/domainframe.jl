# domainframe.jl

"""
A DomainFrame is the restriction of a basis to a subset of its domain. This results
in a frame.
"""
immutable DomainFrame{D,B,N,T} <: AbstractFrame{N,T}
    domain      ::  D
    basis       ::  B

    function DomainFrame(domain::AbstractDomain{N}, basis::FunctionSet{N,T})
        @assert is_basis(basis) == True()
        
        new(domain, basis)
    end
end

DomainFrame{N,T}(domain::AbstractDomain{N}, basis::FunctionSet{N,T}) =
    DomainFrame{typeof(domain),typeof(basis),N,T}(domain, basis)

basis(f::DomainFrame) = f.basis

domain(f::DomainFrame) = f.domain

promote_eltype{S}(f::DomainFrame, ::Type{S}) =
    DomainFrame(f.domain, promote_eltype(f.basis, S))

resize(f::DomainFrame, n) = DomainFrame(domain(f), resize(basis(f), n))


for op in (:size, :length)
    @eval $op(f::DomainFrame, args...) = $op(f.basis, args...)
end

for op in (:differentiation_operator, :antidifferentiation_operator)
    @eval function $op(f::DomainFrame, order)
        op = $op(basis(f), order)
        WrappedOperator(f,DomainFrame(domain(f),dest(op)),op)
    end
end



# Should we check whether x lies in the domain?
call_set(fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion(basis(s), coef, x...)

call_set!(result, fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion!(result, basis(s), coef, x...)


"""
Make a DomainFrame, but match tensor product domains with tensor product sets in a suitable way.

For example: an interval ⊗ a disk (= a cylinder) combined with a 3D Fourier series, leads to a
tensor product of a Fourier series on the interval ⊗ a 2D Fourier series on the disk.
"""
a = 0

function domainframe{TD,DN,LEN}(domain::TensorProductDomain{TD,DN,LEN}, basis::TensorProductSet)
    domainframes = FunctionSet[]
    dc = 1
    for i = 1:LEN
        range = dc:dc+DN[i]-1
        push!(domainframes, DomainFrame(subdomain(domain, i), set(basis, range)))
        dc += DN[i]
    end
    TensorProductSet(domainframes...)
end

domainframe(domain::AbstractDomain, basis::FunctionSet) = DomainFrame(domain, basis)

