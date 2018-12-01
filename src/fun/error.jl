
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

function residual(f::Function, F::DictFun; residualtype = :l2, options...)
    A, B = discretization(dictionary(F), f; options...)
    R = A*coefficients(F)-B
    if residualtype == :l2
        norm(R)
    elseif residualtype == :relativel2
        norm(R)/norm(B)
    elseif residualtype == :maxnorm
        maximum(abs.(R))
    else
        error("Unknown residualtype: $residualtype")
    end
end
