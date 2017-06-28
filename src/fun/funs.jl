# funs.jl

abstract type AbstractFun end

"""
A SetFun corresponds to an expansion in a function set, but it adds a simple user
interface for computing with functions.
"""
struct SetFun{N,T} <: AbstractFun
    expansion   ::  SetExpansion

    SetFun{N,T}(e::SetExpansion) where {N,T} = new(e)
end

SetFun(e::SetExpansion, args...) = SetFun{ndims(e),eltype(e)}(e, args...)

SetFun{N,T}(frame::FunctionSet{N,T}, coefficients = zeros(frame), args...) =
    SetFun{N,T}(SetExpansion(frame, coefficients), args...)

SetFun(domain::Domain, basis::FunctionSet, args...) = SetFun(ExtensionFrame(domain, basis), args...)

SetFun1d{T} = SetFun{1,T}
SetFun2d{T} = SetFun{2,T}
SetFun3d{T} = SetFun{3,T}

expansion(fun::SetFun) = fun.expansion

Base.broadcast(fun::SetFun, x...) = broadcast(expansion(fun), x...)

for op in (:set, :ndims, :coefficients, :eltype, :numtype, :length, :size)
    @eval $op(fun::SetFun) = $op(fun.expansion)
end

for op in (:ctranspose, :∫, :∂x, :∂y, :∂z, :∫∂x, :∫∂y, :∫∂z, :differentiate, :antidifferentiate)
    @eval $op{N,T}(fun::SetFun{N,T}, args...) = SetFun{N,T}($op(fun.expansion, args...))
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

function matrix(fun::SetFun; sampling_factor=2)
    problem = FE_DiscreteProblem(domain(fun), basis(fun), sampling_factor)
    matrix(operator(problem))
end

function sampling_grid(fun::SetFun; sampling_factor=2)
    problem = FE_DiscreteProblem(domain(fun), basis(fun), sampling_factor)
    grid(time_basis_restricted(problem))
end

# Delegate operator applications to the underlying expansion
function (*)(op::AbstractOperator, fun::SetFun)
    @assert src(op) == basis(set(fun))
    SetFun(domain(fun),dest(op),op*coefficients(fun))
end

# Delegate all calling to the underlying expansion.
(fun::SetFun)(x...) = expansion(fun)(x...)


show(io::IO, fun::SetFun) = show(io, fun, set(fun))

function show(io::IO, fun::SetFun, set::FunctionSet)
  println(io, "A ", ndims(fun), "-dimensional SetFun with ", length(coefficients(fun)), " degrees of freedom.")
  println(io, "Basis: ", name(set))
end

function show(io::IO, fun::SetFun, set::ExtensionFrame)
    println(io, "A ", ndims(fun), "-dimensional SetFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(basis(set)))
    println(io, "Domain: ", domain(set))
end

getindex(expansion::SetExpansion, domain::Domain) = restrict(expansion, domain)

getindex(fun::SetFun, domain::Domain) = restrict(expansion(fun), domain)

restrict(expansion::SetExpansion, domain::Domain) = _restrict(expansion, set(expansion), domain)

function _restrict(expansion::SetExpansion, set::ExtensionFrame, domain1::Domain)
    @assert ndims(set) == ndims(domain1)

    domain2 = domain(set)
    newdomain = domain1 ∩ domain2
    SetFun(newdomain, basis(set), coefficients(expansion))
end

function _restrict(expansion::SetExpansion, set::FunctionSet, domain::Domain)
    @assert ndims(set) == ndims(domain)
    # We should check here whether the given domain lies in the support of the set
    SetFun(domain, set, coefficients(expansion))
end

# Get the mean approximation error in random interior points.
function abserror{N}(f::Function,F::SetFun{N};vals::Int=200)
    # Use the bounding box around the domain
    box = boundingbox(domain(F))
    point=Array{numtype(F)}(N)
    elements=0
    error=0
    # Generate some points inside the domain, and compare with the target function
    while elements < vals
        for j in 1:N
            point[j]=left(box)[j]+(right(box)[j]-left(box)[j])*rand(1)[1]
        end
        vpoint = SVector{N}(point)
        if in(vpoint,domain(F))
            elements+=1
            error+=abs(f(vpoint...)-F(vpoint...))
        end
    end
    return error/elements
end

# Get the max approximation error in random interior points
function maxerror{N}(f::Function,F::SetFun{N};vals::Int=200)
    # Use the bounding box around the domain
    box = boundingbox(domain(F))
    point=Array{numtype(F)}(N)
    elements=0
    error=0
    # Generate some points inside the domain, and compare with the target function
    while elements < vals
        for j in 1:N
            point[j]=left(box)[j]+(right(box)[j]-left(box)[j])*rand(1)[1]
        end
        N == 1 ? vpoint = point[1] : vpoint = SVector(point...)
        if in(vpoint,domain(F))
            elements+=1
            error=max(error,abs(f(vpoint...)-F(vpoint...)))
        end
    end
    return error
end

function residual(f::Function, F::SetFun ;  options...)
    op = oversampled_evaluation_operator(basis(F),domain(F); options...)[1]
    rhs = project(dest(op),f)
    # There should be an easier way of getting the inverse of the normalization
    norm(op*coefficients(F)-rhs)
end
