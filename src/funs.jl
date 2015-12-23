# funs.jl

" A FrameFun is an expansion in a basis on a subset of its domain (i.e. in a DomainFrame)."
typealias FrameFun{N,T,D,B,ELT,ID} SetExpansion{DomainFrame{D,B,N,T},ELT,ID}

domain(fun::SetExpansion) = domain(fun, set(fun))
domain(fun::SetExpansion, s::DomainFrame) = domain(s)

basis(fun::FrameFun) = basis(fun, set(fun))
basis(fun::FrameFun, s::DomainFrame) = set(s)


function show_setexpansion(io::IO, fun::SetExpansion, frame::DomainFrame)
    println(io, "A ", dim(fun), "-dimensional FrameFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(basis(frame)))
    println(io, "Domain: ", domain(frame))
end



"""
Construct an FE problem for the given domain, using default values if necessary.
"""
function fe_problem{Basis <: FunctionSet,ELT}(domain, ::Type{Basis}, ::Type{ELT};
    n = default_frame_n(domain, Basis),
    T = default_frame_T(domain, Basis),
    s = default_frame_sampling(domain, Basis),
    solver = default_frame_solver(domain),
    args...)
    
    problem = discretize_problem(domain, n, T, s, Basis, ELT)
    sol = solver(problem)

    (problem, sol)
end

# Detect a suitable element type
fe_problem(domain, Basis) = _fe_problem(domain, Basis, isreal(Basis))
_fe_problem{N,T}(domain::AbstractDomain{N,T}, Basis, isreal::Type{True}) = fe_problem(domain, Basis, T)
_fe_problem{N,T}(domain::AbstractDomain{N,T}, Basis, isreal::Type{False}) = fe_problem(domain, Basis, Complex{T})


ExpFun(f::Function; args...) = Fun(FourierBasis, f; args...)
ExpFun(f::Function, domain; args...) = Fun(FourierBasis, f, domain; args...)

ChebyFun(f::Function; args...) = Fun(ChebyshevBasis, f; args...)
ChebyFun(f::Function, domain; args...) = Fun(ChebyshevBasis, f, domain; args...)


function Fun{Basis <: FunctionSet}(::Type{Basis}, f::Function, domain = default_frame_domain_1d(); args...)
    ELT = eltype(f, domain, Basis)
    (problem,solver) = fe_problem(domain, Basis, ELT; args...)
    solve(solver, f, problem)
end


function eltype{Basis <: AbstractBasis}(f::Function, domain, ::Type{Basis})
    ELT = numtype(domain)
    RT = Base.return_types(f,fill(numtype(domain),dim(domain)))
    if length(RT) > 0
        if isreal(Basis)==Val{false} || (RT[1] <: Complex)
            Complex{ELT}
        else
            ELT
        end
    else
        if isreal(Basis) == Val{false}
            Complex{ELT}
        else
            ELT
        end
    end
end


#call(fun::FrameFun, x...) = in([x[i] for i=1:length(x)], domain(fun)) ? call(fun.expansion, x...) : NaN
call_set(fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion(basis(s), coef, x...)

call_set!(result, fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion!(result, basis(s), coef, x...)


