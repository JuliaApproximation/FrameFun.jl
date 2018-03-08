# funs.jl

abstract type AbstractFun end

"""
A SetFun corresponds to an expansion in a function set, but it adds a simple user
interface for computing with functions.
"""
struct SetFun{T} <: AbstractFun
    expansion   ::  SetExpansion
end

SetFun(e::SetExpansion, args...) = SetFun{domaintype(set(e))}(e, args...)

SetFun{T}(frame::FunctionSet{T}, coefficients = zeros(frame), args...) =
    SetFun{T}(SetExpansion(frame, coefficients), args...)

SetFun(domain::Domain, basis::FunctionSet, args...) = SetFun(ExtensionFrame(domain, basis), args...)

# Warning: not all 2d function sets have SVector{2,T} type, they could have (S,T) type
const SetFun1d{T <: Number} = SetFun{T}
const SetFun2d{T} = SetFun{SVector{2,T}}
const SetFun3d{T} = SetFun{SVector{3,T}}

expansion(fun::SetFun) = fun.expansion

Base.broadcast(fun::SetFun, x...) = broadcast(expansion(fun), x...)

for op in (:set, :dimension, :coefficients, :eltype, :numtype, :length, :size)
    @eval $op(fun::SetFun) = $op(fun.expansion)
end

for op in (:ctranspose, :∫, :∂x, :∂y, :∂z, :∫∂x, :∫∂y, :∫∂z, :differentiate, :antidifferentiate)
    @eval $op{T}(fun::SetFun{T}, args...) = SetFun{T}($op(fun.expansion, args...))
end

for op in (:ExtensionFrame, :basis)
    @eval $op(fun::SetFun) = $op(fun, set(fun))
end

domain(fun::SetFun) = domain(set(fun))

for op in (:+, :-, :*)
    @eval $op(fun1::SetFun,fun2::SetFun) = SetFun($op(fun1.expansion,fun2.expansion))
end

for op in (:+, :-, :*)
    @eval $op(a::Number,fun::SetFun) = SetFun($op(a,fun.expansion))
end

for op in (:+, :-, :*)
    @eval $op(fun::SetFun,a::Number) = $op(a,fun)
end

ExtensionFrame(fun::SetFun, set::ExtensionFrame) = set

domain(fun::SetFun, set::ExtensionFrame) = domain(set)

basis(fun::SetFun, set::ExtensionFrame) = basis(set)

function matrix(fun::SetFun; options...)
    op = oversampled_evaluation_operator(span(basis(fun),eltype(fun)),domain(fun);  options...)[1]
    matrix(op)
end

function sampling_grid(fun::SetFun; sampling_factor=2)
    problem = FE_DiscreteProblem(domain(fun), basis(fun), sampling_factor)
    grid(time_basis_restricted(problem))
end

# Delegate operator applications to the underlying expansion
function (*)(op::AbstractOperator, fun::SetFun)
    @assert src(op) == span(basis(set(fun)))
    SetFun(domain(fun),set(dest(op)),op*coefficients(fun))
end

# Delegate all calling to the underlying expansion.
(fun::SetFun)(x...) = expansion(fun)(x...)


# show(io::IO, fun::SetFun) = show(io, fun, set(fun))

function show(io::IO, fun::SetFun, set::FunctionSet)
  println(io, "A ", dimension(fun), "-dimensional SetFun with ", length(coefficients(fun)), " degrees of freedom.")
  println(io, "Basis: ", name(set))
end

function show(io::IO, fun::SetFun, set::ExtensionFrame)
    println(io, "A ", dimension(fun), "-dimensional SetFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(basis(set)))
    println(io, "Domain: ", domain(set))
end

getindex(expansion::SetExpansion, domain::Domain) = restrict(expansion, domain)

getindex(fun::SetFun, domain::Domain) = restrict(expansion(fun), domain)

restrict(expansion::SetExpansion, domain::Domain) = _restrict(expansion, set(expansion), domain)

function _restrict(expansion::SetExpansion, set::ExtensionFrame, domain1::Domain)
    @assert dimension(set) == dimension(domain1)

    domain2 = domain(set)
    newdomain = domain1 ∩ domain2
    SetFun(newdomain, basis(set), coefficients(expansion))
end

function _restrict(expansion::SetExpansion, set::FunctionSet, domain::Domain)
    @assert dimension(set) == dimension(domain)
    # We should check here whether the given domain lies in the support of the set
    SetFun(domain, set, coefficients(expansion))
end

# Get the mean approximation error in random interior points.
function abserror{N}(f::Function,F::SetFun{N};vals::Int=200)
    rgrid = randomgrid(domain(F),vals)
    Fval = F(rgrid)
    # TODO: type based on space of F
    fval = sample(rgrid,f,eltype(F))
    return sum(abs.(Fval-fval))/vals
end

# Get the max approximation error in random interior points
function maxerror{N}(f::Function,F::SetFun{N};vals::Int=200)
    rgrid = randomgrid(domain(F),vals)
    Fval = F(rgrid)
    # TODO: type based on space of F
    fval = sample(rgrid,f,eltype(F))
    return maximum(abs.(Fval-fval))
end

using QuadGK
function L2error(f::Function, F::SetFun{T}; reltol = eps(real(T)), abstol = 0, options...) where {T}
    I = QuadGK.quadgk(x->abs(F(x)-f(x))^2, left(set(F)), right(set(F)), reltol=reltol, abstol=abstol)
    @assert I[2] < 100max(reltol*I[1],abstol)
    sqrt(I[1])
end

function residual(f::Function, F::SetFun ;  options...)
    op = oversampled_evaluation_operator(span(basis(F)),domain(F); options...)[1]
    rhs = project(dest(op),f)


    norm(op*coefficients(F)-rhs)
end

function relresidual(f::Function, F::SetFun ;  options...)
    op = oversampled_evaluation_operator(span(basis(F)),domain(F); options...)[1]
    rhs = project(dest(op),f)


    norm(op*coefficients(F)-rhs)/norm(rhs)
end

function residualmax(f::Function, F::SetFun ;  options...)
    op = oversampled_evaluation_operator(span(basis(F)),domain(F); options...)[1]
    rhs = project(dest(op),f)


    maximum(abs(op*coefficients(F)-rhs))
end

    

    
