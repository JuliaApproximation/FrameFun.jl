# funs.jl

abstract AbstractFun

"""
A FrameFun corresponds to an expansion in a function set, but it adds a simple user
interface for computing with functions.
"""
immutable FrameFun{N,T} <: AbstractFun
    expansion   ::  SetExpansion
    approx_op   ::  AbstractOperator

    FrameFun(e::SetExpansion, approx_op = approximation_operator(set(e))) = new(e, approx_op)
end

FrameFun(e::SetExpansion, args...) = FrameFun{dim(e),eltype(e)}(e, args...)

FrameFun{N,T}(frame::FunctionSet{N,T}, coefficients = zeros(eltype(frame), size(frame)), args...) =
    FrameFun{N,T}(SetExpansion(frame, coefficients), args...)

FrameFun(domain::AbstractDomain, basis::FunctionSet, args...) = FrameFun(DomainFrame(domain, basis), args...)

expansion(fun::FrameFun) = fun.expansion

for op in (:set, :dim, :coefficients, :eltype, :numtype)
    @eval $op(fun::FrameFun) = $op(fun.expansion)
end

for op in (:ctranspose, :∫, :∂x, :∂y, :∂z, :∫∂x, :∫∂y, :∫∂z, :differentiate, :antidifferentiate)
    @eval $op{N,T}(fun::FrameFun{N,T}, args...) = FrameFun{N,T}($op(fun.expansion, args...))
end

for op in (:domainframe, :domain, :basis)
    @eval $op(fun::FrameFun) = $op(fun, set(fun))
end
    
domainframe(fun::FrameFun, set::DomainFrame) = set

domain(fun::FrameFun, set::DomainFrame) = domain(set)

basis(fun::FrameFun, set::DomainFrame) = basis(set)

function matrix(fun::FrameFun)
    problem = fe_problem(domain(fun), basis(fun), eltype(fun))[1]
    matrix(operator(problem))
end

function sampling_grid(fun::FrameFun)
    problem = fe_problem(domain(fun), basis(fun), eltype(fun))[1]
    grid(time_basis_restricted(problem))
end

# Delegate all calling to the underlying expansion.
call(fun::FrameFun, x...) = call(expansion(fun), x...)


show(io::IO, fun::FrameFun) = show(io, fun, set(fun))

function show(io::IO, fun::FrameFun, set::DomainFrame)
    println(io, "A ", dim(fun), "-dimensional FrameFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(basis(set)))
    println(io, "Domain: ", domain(set))
end

getindex(fun::FrameFun, x...) = getindex(fun, set(fun), x...)

function getindex(fun::FrameFun, set::DomainFrame, domain1::AbstractDomain)
    @assert dim(fun) == dim(domain1)

    domain2 = domain(fun)
    newdomain = domain1 ∩ domain2
    FrameFun(newdomain, basis(fun), coefficients(fun))
end





