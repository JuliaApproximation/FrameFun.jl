# funs.jl

" A FrameFun is an expansion in a basis on a subset of its domain (i.e. in a DomainFrame)."
typealias FrameFun{N,T,D,B,ELT,ID} SetExpansion{DomainFrame{D,B,N,T},ELT,ID}

domain(fun::SetExpansion) = domain(fun, set(fun))
domain(fun::SetExpansion, set::DomainFrame) = domain(set)

basis(fun::FrameFun) = basis(fun, set(fun))
basis(fun::FrameFun, s::DomainFrame) = set(s)

function show{D,B,N}(io::IO, fun::FrameFun{D,B,N})
    println(io, "A ", N, "-dimensional FrameFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(basis(fun)))
    println(io, "Domain: ", domain(fun))
end

function ExpFun(f::Function, domain = default_fourier_domain_1d(),
             solver_type = default_fourier_solver(domain);
             n = default_fourier_n(domain),
             T = default_fourier_T(domain),
             s = default_fourier_sampling(domain))
    Fun(FourierBasis,f,domain,solver_type,n=n,T=T,s=s)
end

function ChebyFun(f::Function, domain = default_fourier_domain_1d(),
             solver_type = default_fourier_solver(domain);
             n = default_fourier_n(domain),
             T = default_fourier_T(domain),
             s = default_fourier_sampling(domain))
    Fun(ChebyshevBasis,f,domain,solver_type,n=n,T=T,s=s)
end

function Fun{Basis <: AbstractBasis}(::Type{Basis}, f::Function, domain = default_fourier_domain_1d(),
             solver_type = default_fourier_solver(domain);
             n = default_fourier_n(domain),
             T = default_fourier_T(domain),
             s = default_fourier_sampling(domain))
    ELT = eltype(f, domain, Basis)
    problem = discretize_problem(domain, n, T, s, Basis, ELT)
    solver = solver_type(problem)
    solve(solver, f, problem, ELT)
end

function eltype{Basis <: AbstractBasis}(f::Function, domain, ::Type{Basis})
    ELT = numtype(domain)
    RT = Base.return_types(f,fill(numtype(domain),dim(domain)))[1]
    if isreal(Basis)==Val{false} || (RT <: Complex)
        Complex{ELT}
    else
        ELT
    end
end


#call(fun::FrameFun, x...) = in([x[i] for i=1:length(x)], domain(fun)) ? call(fun.expansion, x...) : NaN
call_set(fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion(basis(s), coef, x...)

call_set!(result, fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion!(result, basis(s), coef, x...)


