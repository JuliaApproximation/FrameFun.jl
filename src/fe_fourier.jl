# fe_fourier.jl


## "An ExpFun is a FrameFun based on Fourier series."
## ExpFun(f::Function; n=n, T=T, args...) = Fun(FourierBasis, f ; args...)
## # one of three things should be provided: n(tuple), domain or basis
## # only n
## ExpFun(f::Function, n::Int)
## ExpFun{N}(f::Function, n::Ntuple{N}) = Fun(FourierBasis, f, domain; args...)

## "A ChebyFun is a FrameFun based on Chebyshev polynomials."
## ChebyFun(f::Function; args...) = Fun(ChebyshevBasis, f; args...)
## ChebyFun(f::Function, domain; args...) = Fun(ChebyshevBasis, f, domain; args...)


function Fun(f::Function, Basis::FunctionSet, domain = default_frame_domain_1d(typeof(Basis)); args...)
    ELT = eltype(f, Basis)
    (problem,solver) = fe_problem(domain, Basis, ELT; args...)
    coef = solve(solver, f, problem)
    FrameFun(domain, frequency_basis(problem), coef)
end


function eltype(f::Function, Basis)
    ELT = numtype(Basis)
    RT = Base.return_types(f,fill(numtype(Basis),dim(Basis)))
    if length(RT) > 0
        if isreal(typeof(Basis)) == False || (RT[1] <: Complex)
            Complex{ELT}
        else
            ELT
        end
    else
        if isreal(typeof(Basis)) == False
            Complex{ELT}
        else
            ELT
        end
    end
end



## """
## Construct an FE problem for the given domain, using default values if necessary.
## """
## remove
function fe_problem(domain, Basis, ELT;
    s = 2,
    solver = default_frame_solver(domain, Basis),
    args...)
    problem = discretize_problem(domain, Basis, ELT, s)
    sol = solver(problem; args...)

    (problem, sol)
end




function discretize_problem{N,T}(domain::AbstractDomain{N,T}, fbasis1, ELT, st)

    rgrid, fbasis2 = oversampled_grid(domain, fbasis1, st)
    grid1 = grid(fbasis1)
    grid2 = grid(fbasis2)

    tbasis1 = DiscreteGridSpace(grid1, ELT)
    tbasis2 = DiscreteGridSpace(grid2, ELT)
    tbasis_restricted = DiscreteGridSpace(rgrid, ELT)

    FE_DiscreteProblem(domain, fbasis1, fbasis2, tbasis1, tbasis2, tbasis_restricted)
end

function discretize_problem{TD,DN,LEN,N,T}(domain::TensorProductDomain{TD,DN,LEN,N,T}, Basis, ELT, st)
    problems = FE_Problem[]
    dc = 1
    for i = 1:LEN
        range = dc:dc+DN[i]-1
        p = discretize_problem(subdomain(domain,i), set(Basis,range), ELT, st)
        push!(problems, p)
        dc += DN[i]
    end
    problem = FE_TensorProductProblem(problems...)
end


######################
# Default parameters
######################

# The default parameters are:
# - n: number of degrees of freedom in the approximation
# - T: the extension parameter (size of the extended domain vs the original domain)
# - sampling: the (over)sampling factor
# - domain_nd: the domain
# - solver: the solver to use for solving the FE problem


# The default basis is a FourierBasis
for op in (:default_frame_1d, :default_frame_2d, :default_frame_3d)
    @eval $op() = $op(FourierBasis)
end

for op in (:default_frame_n, :default_frame_T, :default_frame_sampling, :default_frame_solver)
    @eval $op(domain) = $op(domain, FourierBasis)
end

default_frame_domain_1d(Basis) = Interval()
default_frame_domain_2d(Basis) = Interval()⊗Interval()
default_frame_domain_3d(Basis) = Interval() ⊗ Interval() ⊗ Interval()


default_frame_n(domain::AbstractDomain1d, Basis) = 61
default_frame_n(domain::AbstractDomain2d, Basis) = (31, 31)
default_frame_n(domain::AbstractDomain3d, Basis) = (7, 7, 7)


function default_frame_n(domain::TensorProductDomain, Basis)
    s = [default_frame_n(domainlist(domain)[1], Basis)...]
    for i = 2:tp_length(domain)
        s = [s; default_frame_n(domainlist(domain)[i], Basis)...]
    end
    s = round(Int,s/dim(domain))
    tuple(s...)
end


default_frame_T{T}(domain::AbstractDomain{1,T}, Basis) = 2*one(T)
default_frame_T{N,T}(domain::AbstractDomain{N,T}, Basis) = ntuple(i->2*one(T),N)

default_frame_solver(domain, Basis) = FE_ProjectionSolver

default_frame_solver{N}(domain::AbstractDomain, Basis::FunctionSet{N,BigFloat}) = FE_DirectSolver
default_frame_solver{N}(domain::AbstractDomain, Basis::FunctionSet{N,Complex{BigFloat}}) = FE_DirectSolver

