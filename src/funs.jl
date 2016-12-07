# funs.jl

abstract AbstractFun

"""
A FrameFun corresponds to an expansion in a function set, but it adds a simple user
interface for computing with functions.
"""
immutable FrameFun{N,T} <: AbstractFun
    expansion   ::  SetExpansion

    FrameFun(e::SetExpansion) = new(e)
end

FrameFun(e::SetExpansion, args...) = FrameFun{ndims(e),eltype(e)}(e, args...)

FrameFun{N,T}(frame::FunctionSet{N,T}, coefficients = zeros(eltype(frame), size(frame)), args...) =
    FrameFun{N,T}(SetExpansion(frame, coefficients), args...)

FrameFun(domain::AbstractDomain, basis::FunctionSet, args...) = FrameFun(DomainFrame(domain, basis), args...)

typealias FrameFun1d{T} FrameFun{1,T}
typealias FrameFun2d{T} FrameFun{2,T}
typealias FrameFun3d{T} FrameFun{3,T}

expansion(fun::FrameFun) = fun.expansion

for op in (:set, :ndims, :coefficients, :eltype, :numtype)
    @eval $op(fun::FrameFun) = $op(fun.expansion)
end

for op in (:ctranspose, :∫, :∂x, :∂y, :∂z, :∫∂x, :∫∂y, :∫∂z, :differentiate, :antidifferentiate)
    @eval $op{N,T}(fun::FrameFun{N,T}, args...) = FrameFun{N,T}($op(fun.expansion, args...))
end

for op in (:domainframe, :domain, :basis)
    @eval $op(fun::FrameFun) = $op(fun, set(fun))
end

for op in (:+, :-, :*)
    @eval $op(fun1::FrameFun,fun2::FrameFun) = FrameFun($op(fun1.expansion,fun2.expansion))
end

for op in (:+, :-, :*)
    @eval $op(a::Number,fun::FrameFun) = FrameFun($op(a,fun.expansion))
end

for op in (:+, :-, :*)
    @eval $op(fun::FrameFun,a::Number) = $op(a,fun)
end

domainframe(fun::FrameFun, set::DomainFrame) = set

domain(fun::FrameFun, set::DomainFrame) = domain(set)

basis(fun::FrameFun, set::DomainFrame) = basis(set)

function matrix(fun::FrameFun; sampling_factor=2)
    problem = FE_DiscreteProblem(domain(fun), basis(fun), sampling_factor)
    matrix(operator(problem))
end

function sampling_grid(fun::FrameFun; sampling_factor=2)
    problem = FE_DiscreteProblem(domain(fun), basis(fun), sampling_factor)
    grid(time_basis_restricted(problem))
end

# Delegate operator applications to the underlying expansion
function (*)(op::AbstractOperator, fun::FrameFun)
    @assert src(op) == basis(set(fun))
    FrameFun(domain(fun),dest(op),op*coefficients(fun))
end

# Delegate all calling to the underlying expansion.
@compat (fun::FrameFun)(x...) = expansion(fun)(x...)


show(io::IO, fun::FrameFun) = show(io, fun, set(fun))

function show(io::IO, fun::FrameFun, set::FunctionSet)
  println(io, "A ", ndims(fun), "-dimensional FrameFun with ", length(coefficients(fun)), " degrees of freedom.")
  println(io, "Basis: ", name(set))
end

function show(io::IO, fun::FrameFun, set::DomainFrame)
    println(io, "A ", ndims(fun), "-dimensional FrameFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(basis(set)))
    println(io, "Domain: ", domain(set))
end

getindex(fun::FrameFun, x...) = getindex(fun, set(fun), x...)

function getindex(fun::FrameFun, set::DomainFrame, domain1::AbstractDomain)
    @assert ndims(fun) == ndims(domain1)

    domain2 = domain(fun)
    newdomain = domain1 ∩ domain2
    FrameFun(newdomain, basis(fun), coefficients(fun))
end

# Get the mean approximation error in random interior points.
function abserror{N}(f::Function,F::FrameFun{N};vals::Int=200)
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
        N == 1 ? vpoint = point[1] : vpoint = Vec(point...)
        if in(vpoint,domain(F))
            elements+=1
            error+=abs(f(vpoint...)-F(vpoint...))
        end
    end
    return error/elements
end

# Get the max approximation error in random interior points
function maxerror{N}(f::Function,F::FrameFun{N};vals::Int=200)
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
        N == 1 ? vpoint = point[1] : vpoint = Vec(point...)
        if in(vpoint,domain(F))
            elements+=1
            error=max(error,abs(f(vpoint...)-F(vpoint...)))
        end
    end
    return error
end

function residual(f::Function, F::FrameFun ; sampling_factor=2, options...)
    problem = FE_DiscreteProblem(domain(F), basis(F), sampling_factor; options...)
    op = operator(problem)
    rhs = sample(grid(dest(op)), f, eltype(src(op)))
    op = full_transform_operator(src(op),dest(op))
    norm(op*coefficients(F)-rhs)
end
