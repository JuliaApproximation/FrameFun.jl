# domainframe.jl

# A DomainFrame is the restriction of a basis to a subset of its domain. This results
# in a frame.
immutable DomainFrame{D,B,N,T} <: AbstractFrame{N,T}
    domain      ::  D
    basis       ::  B
end

DomainFrame{N,T}(domain::AbstractDomain{N,T},basis::AbstractBasis{N,T}) = 
    DomainFrame{typeof(domain),typeof(basis),N,T}(domain,basis)

basis(f::DomainFrame) = f.basis

domain(f::DomainFrame) = f.domain


