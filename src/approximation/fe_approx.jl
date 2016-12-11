# fe_fourier.jl


## "An ExpFun is a SetFun based on Fourier series."
## ExpFun(f::Function; n=n, T=T, args...) = Fun(FourierBasis, f ; args...)
## # one of three things should be provided: n(tuple), domain or basis
## # only n

## ExpFun(f::Function, n::Int)
## ExpFun{N}(f::Function, n::Ntuple{N}) = Fun(FourierBasis, f, domain; args...)

## "A ChebyFun is a SetFun based on Chebyshev polynomials."
## ChebyFun(f::Function; args...) = Fun(ChebyshevBasis, f; args...)
## ChebyFun(f::Function, domain; args...) = Fun(ChebyshevBasis, f, domain; args...)


function Fun(f::Function, basis::FunctionSet, domain::AbstractDomain; options...)
    ELT = eltype(f, basis)
    frame = domainframe(domain, promote_eltype(basis, ELT))
    A = approximation_operator(frame; options...)
    coef = A * f
    SetFun(domain, dest(A), coef)
end

# We assume f as a function is type stable.
function eltype(f::Function, basis)
    ELT = eltype(basis)
    # We only test for the return type in zero
    RT = typeof(f(fill(zero(ELT),ndims(basis))...))
    if (RT <: Complex)
        complexify(ELT)
    else
        ELT
    end
end



function approximation_operator(set::DomainFrame;
    sampling_factor = 2, solver = default_frame_solver(domain(set), basis(set)), options... )

    problem = FE_DiscreteProblem(domain(set), basis(set), sampling_factor; options...)
    solver(problem; options...)
end


immutable FE_BestSolver
end

function FE_BestSolver(problem::FE_DiscreteProblem; options...)
    if has_transform(frequency_basis(problem))
        R = estimate_plunge_rank(problem)
        if R < size(problem, 2)/2
            FE_ProjectionSolver(problem; options...)
        else
            FE_DirectSolver(problem; options...)
        end
    else
        # Don't bother with a fast algorithm if there is no fast transform
        FE_DirectSolver(problem; options...)
    end
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


default_frame_domain_1d(basis) = Interval()
default_frame_domain_2d(basis) = Interval() ⊗ Interval()
default_frame_domain_3d(basis) = Interval() ⊗ Interval() ⊗ Interval()


default_frame_n(domain::AbstractDomain1d, basis) = 61
default_frame_n(domain::AbstractDomain2d, basis) = (31, 31)
default_frame_n(domain::AbstractDomain3d, basis) = (7, 7, 7)


function default_frame_n(domain::TensorProductDomain, basis)
    s = [default_frame_n(element(domain,1), basis)...]
    for i = 2:composite_length(domain)
        s = [s; default_frame_n(element(domain,i), basis)...]
    end
    s = round(Int,s/ndims(domain))
    tuple(s...)
end


default_frame_T(domain::AbstractDomain{1}, basis) = 2.0
default_frame_T{N}(domain::AbstractDomain{N}, basis) = ntuple(i->2.0,N)

default_frame_solver(domain, basis) = FE_BestSolver

default_frame_solver{N}(domain::AbstractDomain, basis::FunctionSet{N,BigFloat}) = FE_DirectSolver
default_frame_solver{N}(domain::AbstractDomain, basis::FunctionSet{N,Complex{BigFloat}}) = FE_DirectSolver
