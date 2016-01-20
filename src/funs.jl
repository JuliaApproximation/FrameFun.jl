# funs.jl

abstract AbstractFun

"""
A FrameFun corresponds to an expansion in a function set, but it adds a simple user
interface for computing with functions.
"""
immutable FrameFun{N,T} <: AbstractFun
    expansion   ::  SetExpansion
end


FrameFun{N,T,ELT}(domain::AbstractDomain{N,T}, basis::FunctionSet{N,ELT}, coefficients) =
    FrameFun{N,ELT}(SetExpansion(DomainFrame(domain, basis), coefficients))

FrameFun(domain::AbstractDomain, basis::FunctionSet) = FrameFun(domain, basis, zeros(eltype(basis), size(basis)))

expansion(fun::FrameFun) = fun.expansion

for op in (:set, :dim, :coefficients, :eltype, :numtype)
    @eval $op(fun::FrameFun) = $op(fun.expansion)
end

for op in (:domainframe, :domain, :basis)
    @eval $op(fun::FrameFun) = $op(fun, set(fun))
end

domainframe(fun::FrameFun, set::DomainFrame) = set

domain(fun::FrameFun, set::DomainFrame) = domain(set)

basis(fun::FrameFun, set::DomainFrame) = basis(set)


# Delegate all calling to the underlying expansion.
call(fun::FrameFun, x...) = call(expansion(fun), x...)


show(io::IO, fun::FrameFun) = show(io, fun, set(fun))

@debug function show(io::IO, fun::FrameFun, set::DomainFrame)
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





