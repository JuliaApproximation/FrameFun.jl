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

for op in (:size, :length, :differentiate)
    @eval $op(f::DomainFrame, args...) = $op(f.basis, args...)
end


# Should we check whether x lies in the domain?
call_set(fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion(basis(s), coef, x...)

call_set!(result, fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion!(result, basis(s), coef, x...)



"""
Compute a grid of a larger basis, but restricted to the given domain, using oversampling by the given factor
(approximately) in each dimension.
The result is the tuple (oversampled_grid, larger_basis)
"""
oversampled_grid(set::DomainFrame, sampling_factor = 2) = compute_subgrid(domain(set), basis(set), sampling_factor)

function oversampled_grid(domain::AbstractDomain, basis::FunctionSet, sampling_factor)
    N = dim(basis)
    n_goal = length(basis) * sampling_factor^N
    grid1 = grid(basis)
    grid2 = subgrid(grid1, domain)
    ratio = length(grid2) / length(grid1)
    # This could be way off if the original size was small.
    n = approx_length(basis, ceil(Int, n_goal/ratio))
    large_basis = similar(basis, n)
    grid3 = grid(large_basis)
    grid4 = subgrid(grid3, domain)
    grid4, large_basis
end


