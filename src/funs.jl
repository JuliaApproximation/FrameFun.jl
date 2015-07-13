# funs.jl


###############################################################################################
# Some common properties of Function objects are grouped in AbstractFun
###############################################################################################
abstract AbstractFun{N,T <: FloatingPoint}

dim{N,T}(::AbstractFun{N,T}) = N
dim{N,T}(::Type{AbstractFun{N,T}}) = N
dim{F <: AbstractFun}(::Type{F}) = dim(super(F))

numtype{N,T}(::AbstractFun{N,T}) = T
numtype{N,T}(::Type{AbstractFun{N,T}}) = T
numtype{F <: AbstractFun}(::Type{F}) = numtype(super(F))


# A Fun groups a domain and an expansion
immutable Fun{N,T} <: AbstractFun{N,T}
    domain      ::  AbstractDomain{N,T}
    expansion   ::  SetExpansion        # the frame coefficients
#    problem     ::  FE_Problem{N,T}     # properties of the Fourier extension problem, store for convenience

    function Fun(domain, expansion)
        @assert numtype(expansion) == numtype(domain)
        @assert dim(expansion) == dim(domain)
        
        new(domain, expansion)
    end
end
# Obsolete constructor
## Fun{N,T}(domain::AbstractDomain{N,T}, expansion, problem) = Fun{N,T}(domain, expansion, problem)

Fun{N,T}(domain::AbstractDomain{N,T}, expansion) = Fun{N,T}(domain, expansion)

domain(fun::Fun) = fun.domain

expansion(fun::Fun) = fun.expansion

set(fun::Fun) = set(expansion(fun))

coefficients(fun::Fun) = coefficients(expansion(fun))

#problem(fun::Fun) = fun.problem

function show{N}(io::IO, fun::Fun{N})
    println(io, "A ", N, "-dimensional FrameFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(set(expansion(fun))))
    println(io, "Domain: ", domain(fun))
end

function ExpFun(f::Function, domain = default_fourier_domain_1d(),
        solver_type = default_fourier_solver(domain);
        n = default_fourier_n(domain), T = default_fourier_T(domain),
        s = default_fourier_sampling(domain))

    problem = default_fourier_problem(domain, n, T, s)
    solver = solver_type(problem)

    expansion = solve(solver, f)
    Fun(domain, expansion)
end

# Funs with one solver_type take that as the default
function ExpFun{TD,DN,ID,N}(f::Function, domain::TensorProductDomain{TD,DN,ID,N},
        solver_type = default_fourier_solver(domain);
        n = default_fourier_n(domain), T = default_fourier_T(domain),
        s = default_fourier_sampling(domain))
    problems=FE_DiscreteProblem[]
    dc=1
    for i=1:ID
        push!(problems,default_fourier_problem(subdomain(domain,i),n[dc:dc+DN[i]-1],T[dc:dc+DN[i]-1],s[dc:dc+DN[i]-1]))
        dc=dc+DN[i]
    end
    solver = FE_TensorProductSolver(problems,ntuple(i->solver_type,ID))
    expansion = solve(solver, f)
    Fun(domain, expansion)
end

function ExpFun{TD,DN,ID,N}(f::Function, domain::TensorProductDomain{TD,DN,ID,N},
        solver_types::Tuple;
        n = default_fourier_n(domain), T = default_fourier_T(domain),
        s = default_fourier_sampling(domain))
    problems=FE_DiscreteProblem[]
    dc=1
    for i=1:ID
        push!(problems,default_fourier_problem(subdomain(domain,i),n[dc:dc+DN[i]-1],T[dc:dc+DN[i]-1],s[dc:dc+DN[i]-1]))
        dc=dc+DN[i]
    end
    solver = FE_TensorProductSolver(problems,solver_types)
    expansion = solve(solver, f)
    Fun(domain, expansion)
end


call(fun::Fun, x...) = in([x[i] for i=1:length(x)], domain(fun)) ? call(fun.expansion, x...) : NaN





