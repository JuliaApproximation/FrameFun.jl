# funs.jl

abstract type AbstractFun
end

"""
A `DictFun` corresponds to an expansion in a dictionary, but it adds a simple user
interface for computing with functions.
"""
struct DictFun{S,T} <: AbstractFun
    expansion   ::  Expansion
end

DictFun(e::Expansion, args...) = DictFun{domaintype(dictionary(e)),codomaintype(dictionary(e))}(e, args...)

DictFun(frame::Dictionary{S,T}, coefficients = zeros(frame), args...) where {S,T} =
    DictFun{S,T}(Expansion(frame, coefficients), args...)

DictFun(domain::Domain, basis::Dictionary, args...) = DictFun(ExtensionFrame(domain, basis), args...)

# Warning: not all 2d function sets have SVector{2,T} type, they could have (S,T) type
const DictFun1d{S <: Number,T} = DictFun{S,T}
const DictFun2d{S <: Number,T} = DictFun{SVector{2,S},T}
const DictFun3d{S <: Number,T} = DictFun{SVector{3,S},T}
const DictFun4d{S <: Number,T} = DictFun{SVector{4,S},T}

expansion(fun::DictFun) = fun.expansion

Base.broadcast(fun::DictFun, x...) = broadcast(expansion(fun), x...)

for op in (:dictionary, :dimension, :coefficients, :eltype, :numtype, :length, :size)
    @eval $op(fun::DictFun) = $op(fun.expansion)
end

for op in (:adjoint, :∫, :∂x, :∂y, :∂z, :∫∂x, :∫∂y, :∫∂z, :differentiate, :antidifferentiate)
    @eval $op(fun::DictFun{S,T}, args...) where {S,T} = DictFun{S,T}($op(fun.expansion, args...))
end

for op in (:ExtensionFrame, :basis, :basisspan)
    @eval $op(fun::DictFun) = $op(fun, dictionary(fun))
end

domain(fun::DictFun) = domain(fun, dictionary(fun))

for op in (:+, :-, :*)
    @eval $op(fun1::DictFun,fun2::DictFun) = DictFun($op(fun1.expansion,fun2.expansion))
end

for op in (:+, :-, :*)
    @eval $op(a::Number,fun::DictFun) = DictFun($op(a,fun.expansion))
end

for op in (:+, :-, :*)
    @eval $op(fun::DictFun,a::Number) = $op(a,fun)
end

ExtensionFrame(fun::DictFun, set::ExtensionFrame) = set

domain(fun::DictFun, set::ExtensionFrame) = domain(set)

basis(fun::DictFun, set::ExtensionFrame) = basis(set)

basis(fun::DictFun, dict::Dictionary) = dict

domain(fun::DictFun, dict::Dictionary) = support(dict)

function matrix(fun::DictFun; options...)
    op = oversampled_evaluation_operator(basis(fun),domain(fun);  options...)[1]
    matrix(op)
end

# Delegate operator applications to the underlying expansion
function (*)(op::DictionaryOperator, fun::DictFun)
    @assert (src(op) == dictionary(fun)) | (src(op) == basis(dictionary(fun)))
    DictFun(dest(op),op*coefficients(fun))
end

# Delegate all calling to the underlying expansion.
(fun::DictFun)(x...) = expansion(fun)(x...)


# show(io::IO, fun::DictFun) = show(io, fun, set(fun))

function show(io::IO, fun::DictFun, set::Dictionary)
    println(io, "A ", dimension(fun), "-dimensional DictFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(set))
end

function show(io::IO, fun::DictFun, set::ExtensionFrame)
    println(io, "A ", dimension(fun), "-dimensional DictFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(basis(set)))
    println(io, "Domain: ", domain(set))
end

getindex(expansion::Expansion, domain::Domain) = restrict(expansion, domain)

getindex(fun::DictFun, domain::Domain) = restrict(expansion(fun), domain)

restrict(expansion::Expansion, domain::Domain) = _restrict(expansion, dictionary(expansion), domain)

function _restrict(expansion::Expansion, set::ExtensionFrame, domain1::Domain)
    @assert dimension(set) == dimension(domain1)

    domain2 = domain(set)
    newdomain = domain1 ∩ domain2
    DictFun(newdomain, basis(set), coefficients(expansion))
end

function _restrict(expansion::Expansion, set::Dictionary, domain::Domain)
    @assert dimension(set) == dimension(domain)
    # We should check here whether the given domain lies in the support of the set
    DictFun(domain, set, coefficients(expansion))
end

# Get the mean approximation error in random interior points.
function abserror(f::Function,F::DictFun;vals::Int=200)
    rgrid = randomgrid(domain(F),vals)
    Fval = F(rgrid)
    # TODO: type based on space of F
    fval = sample(rgrid,f,eltype(F))
    return sum(abs.(Fval-fval))/vals
end

# Get the max approximation error in random interior points
function maxerror(f::Function,F::DictFun;vals::Int=200)
    rgrid = randomgrid(domain(F),vals)
    Fval = F(rgrid)
    # TODO: type based on space of F
    fval = sample(rgrid,f,eltype(F))
    return maximum(abs.(Fval-fval))
end

using QuadGK

function L2error(f::Function, F::DictFun{S,T}; rtol = eps(real(T)), atol = 0, options...) where {S,T}
    I = QuadGK.quadgk(x->abs(F(x)-f(x))^2, left(dictionary(F)), right(dictionary(F)), rtol=rtol, atol=atol)
    @assert I[2] < 100max(rtol*I[1],atol)
    sqrt(I[1])
end

function residual(f::Function, F::DictFun ;  options...)
    op = oversampled_evaluation_operator(basis(F),domain(F); options...)[1]
    rhs = project(dest(op),f)
    norm(op*coefficients(F)-rhs)
end

function relresidual(f::Function, F::DictFun ;  options...)
    op = oversampled_evaluation_operator(basis(F),domain(F); options...)[1]
    rhs = project(dest(op),f)


    norm(op*coefficients(F)-rhs)/norm(rhs)
end

function residualmax(f::Function, F::DictFun ;  options...)
    op = oversampled_evaluation_operator(basis(F),domain(F); options...)[1]
    rhs = project(dest(op),f)


    maximum(abs(op*coefficients(F)-rhs))
end
