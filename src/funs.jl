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

immutable Fun{N,T,D,S,ELT,ID} <: AbstractFun{N,T}
    domain      ::  D
    expansion   ::  SetExpansion{S,ELT,ID}     # the frame coefficients

    Fun(domain::AbstractDomain{N,T}, coef::SetExpansion{S,ELT,ID}) = new(domain,coef)
end

Fun{N,T,S,ELT,ID}(domain::AbstractDomain{N,T}, expansion::SetExpansion{S,ELT,ID}) = Fun{N,T,typeof(domain),S,ELT,ID}(domain, expansion)

domain(fun::Fun) = fun.domain

expansion(fun::Fun) = fun.expansion

function ExpFun(f::Function, domain = default_fourier_domain_1d(), problem = default_fourier_problem_1d(), solver = default_fourier_solver(problem))
    expansion = solve(solver, f)
    Fun(domain, expansion)
end


call(fun::Fun, x) = in(x, domain(fun)) ? call(fun.expansion, x...) : NaN


