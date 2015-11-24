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

function Fun(Basis::DataType, f::Function, domain = default_fourier_domain_1d(),
             solver_type = default_fourier_solver(domain);
             n = default_fourier_n(domain),
             T = default_fourier_T(domain),
             s = default_fourier_sampling(domain))
    ELT=Base.return_types(f,fill(numtype(domain),dim(domain)))[1]
    if isreal(Basis)==Val{false}() && (T<:Real)
        ELT=Complex{T}
    end        
    problem = discretize_problem(domain, n, T, s, Basis, ELT)
    solver = solver_type(problem)

    solve(solver, f)
end

function Fun{TD,DN,ID,N}(Basis::DataType, f::Function, domain::TensorProductDomain{TD,DN,ID,N},
        solver_types::Tuple;
        n = default_fourier_n(domain), T = default_fourier_T(domain),
        s = default_fourier_sampling(domain))
    problems=FE_DiscreteProblem[]
    dc=1
    ELT=Base.return_types(f,fill(numtype(domain),dim(domain)))[1]
    if isreal(Basis)==Val{false}() && (T<:Real)
        T=Complex{T}
    end        
    for i=1:ID
        push!(problems,discretize_problem(subdomain(domain,i),n[dc:dc+DN[i]-1],T[dc:dc+DN[i]-1],s[dc:dc+DN[i]-1],Basis,ELT))
        dc=dc+DN[i]
    end
    solver = FE_TensorProductSolver(problems,solver_types)
    solve(solver, f)
end

# Funs with one solver_type take that as the default
function Fun{TD,DN,ID,N}(Basis::DataType, f::Function, domain::TensorProductDomain{TD,DN,ID,N},
        solver_type = default_fourier_solver(domain);
        n = default_fourier_n(domain), T = default_fourier_T(domain),
        s = default_fourier_sampling(domain))
    problems=FE_DiscreteProblem[]
    dc=1
    ELT=Base.return_types(f,fill(numtype(domain),dim(domain)))[1]
    if isreal(Basis)==Val{false}() && (T<:Real)
        ELT=Complex{ELT}
    end        
    for i=1:ID
        push!(problems,discretize_problem(subdomain(domain,i),n[dc:dc+DN[i]-1],T[dc:dc+DN[i]-1],s[dc:dc+DN[i]-1],Basis,ELT))
        dc=dc+DN[i]
    end
    solver = FE_TensorProductSolver(problems,ntuple(i->solver_type,ID))
    solve(solver, f)
end


#call(fun::FrameFun, x...) = in([x[i] for i=1:length(x)], domain(fun)) ? call(fun.expansion, x...) : NaN
call_set(fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion(basis(s), coef, x...)

call_set!(result, fun::SetExpansion, s::DomainFrame, coef, x...) = call_expansion!(result, basis(s), coef, x...)


