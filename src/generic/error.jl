
# Get the mean approximation error in random interior points.
function abserror(f::Function, F::Expansion; vals = 200)
    rgrid = randomgrid(support(F),vals)
    @debug "Random evaluation in abserror"
    Fval = F(rgrid)
    fval = sample(rgrid,f,codomaintype(F))
    sum(abs.(Fval-fval))/vals
end

# Get the max approximation error in random interior points
function maxerror(f::Function, F::Expansion; vals = 200, options...)
    rgrid = randomgrid(support(F),vals)
    @debug "Random evaluation in maxerror"
    Fval = F(rgrid; options...)
    fval = sample(rgrid,f,codomaintype(F))
    maximum(abs.(Fval-fval))
end

function L2error(f::Function, F::Expansion{S,T}; rtol = eps(real(T)), atol = eps(real(T)), options...) where {S,T}
    I = QuadGK.quadgk(x->abs(F(x)-f(x))^2, infimum(support(dictionary(F))), supremum(support(dictionary(F))), rtol=rtol, atol=atol)
    @assert I[2] < 100max(rtol*I[1],atol)
    sqrt(I[1])
end

function residual(f::Function, F::Expansion; residualtype = :l2, options...)
    A, B = full_discretization(f, dictionary(F); options...)
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
