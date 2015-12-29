# domainframe.jl

"""
A DomainFrame is the restriction of a basis to a subset of its domain. This results
in a frame.
"""
immutable DomainFrame{D,B,N,T} <: AbstractFrame{N,T}
    domain      ::  D
    basis       ::  B

    function DomainFrame(domain::AbstractDomain{N,T}, basis::FunctionSet{N,T})
        @assert is_basis(basis) == True()
        
        new(domain, basis)
    end
end

DomainFrame{N,T}(domain::AbstractDomain{N,T}, basis::FunctionSet{N,T}) = 
    DomainFrame{typeof(domain),typeof(basis),N,T}(domain, basis)

basis(f::DomainFrame) = f.basis

domain(f::DomainFrame) = f.domain

for op in (:size, :length)
    @eval $op(f::DomainFrame, args...) = $op(f.basis, args...)
end

# Should we check whether x lies in the domain?
call_set(fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion(basis(s), coef, x...)

call_set!(result, fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion!(result, basis(s), coef, x...)


