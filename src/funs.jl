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
immutable Fun{N,T,D,S,ELT,ID} <: AbstractFun{N,T}
    domain      ::  D
    expansion   ::  SetExpansion{S,ELT,ID}     # the frame coefficients

    Fun(domain::AbstractDomain{N,T}, expansion::SetExpansion{S,ELT,ID}) = new(domain, expansion)
end

Fun{N,T,S,ELT,ID}(domain::AbstractDomain{N,T}, expansion::SetExpansion{S,ELT,ID}) = Fun{N,T,typeof(domain),S,ELT,ID}(domain, expansion)

domain(fun::Fun) = fun.domain

expansion(fun::Fun) = fun.expansion

set(fun::Fun) = set(expansion(fun))

coefficients(fun::Fun) = coefficients(expansion(fun))

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


call(fun::Fun, x...) = in([x[i] for i=1:length(x)], domain(fun)) ? call(fun.expansion, x...) : NaN





