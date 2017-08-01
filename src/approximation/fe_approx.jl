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


function Fun(f::Function, basis::FunctionSet, domain::Domain; options...)
    ELT = real(rangetype(f, basis))
    frame = extensionframe(domain, promote_domaintype(basis, ELT))
    A = approximation_operator(span(frame); options...)
    coef = A * f
    SetFun(domain, set(dest(A)), coef)
end

function fe_problem(basis, domain, sampling_factor; options...)
    frame = ExtensionFrame(domain, basis)
    FE_DiscreteProblem(domain, basis, sampling_factor; options...)
end

function fe_solver(basis, domain; options...)
    frame = ExtensionFrame(domain, basis)
    approximation_operator(span(frame); options...)
end

# We assume f as a function is type stable.
function rangetype(f::Function, basis)
    ELT = rangetype(basis)
    # We only test for the return type in zero
    RT = typeof(f(fill(zero(domaintype(basis)),dimension(basis))...))
    if (RT <: Complex)
        complex(ELT)
    else
        ELT
    end
end

function oversampled_evaluation_operator(S::Span, D::Domain; sampling_factor=2, incboundary=false, options...)
    B = primaryspan(S)
    # Establish time domain grid
    G, lB = oversampled_grid(D,set(B),sampling_factor)

    op = grid_evaluation_operator(S,gridspace(B,G),G)
    # Add boundary points if necessary
    if incboundary
        BG = boundary(grid(lB), D)
        op = [op; grid_evaluation_operator(S,gridspace(B,BG),BG)]
    end
    (op,length(lB))
end

function discrete_approximation_operator(set::ExtensionSpan; solver = default_frame_solver(domain(set), basisspan(set)), options...)
    (op, scaling) = oversampled_evaluation_operator(basisspan(set),domain(set);options...)
    solver(op, scaling; options...)
end

primaryspan(span::Span) = span
function primaryspan(span::BasisFunctions.MultiSetSpan)
    elements(span)[findmax(map(length,elements(span)))[2]]
end

struct FE_BestSolver
end

function FE_BestSolver(op::AbstractOperator, scaling; verbose= false, options...)
    if has_transform(src(op))
        R = estimate_plunge_rank(op)
        if R < size(op, 2)/2
            verbose && println("Estimated plunge rank $R smaller than $(size(op,2))/2 -> projection solver ")
            FE_ProjectionSolver(op, scaling; verbose=verbose,options...)
        else
            verbose && println("Estimated plunge rank $R greater than $(size(op,2))/2 -> direct solver ")
            FE_DirectSolver(op, scaling; verbose=verbose,options...)
        end
    else
        # Don't bother with a fast algorithm if there is no fast transform
        FE_DirectSolver(op, scaling; verbose=verbose, options...)
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


default_frame_n(domain::Domain1d, basis) = 61
default_frame_n(domain::Domain2d, basis) = (31, 31)
default_frame_n(domain::Domain3d, basis) = (7, 7, 7)


function default_frame_n(domain::ProductDomain, basis)
    s = [default_frame_n(element(domain,1), basis)...]
    for i = 2:nb_elements(domain)
        s = [s; default_frame_n(element(domain,i), basis)...]
    end
    s = round.(Int,s/dimension(domain))
    tuple(s...)
end


# default_frame_T(domain::Domain{1}, basis) = 2.0
# default_frame_T{N}(domain::Domain{N}, basis) = ntuple(i->2.0,N)

default_frame_solver(domain, basis) = FE_BestSolver

default_frame_solver(domain::Domain, basis::FunctionSet{SVector{N,BigFloat}}) where {N} = FE_DirectSolver
default_frame_solver(domain::Domain, basis::FunctionSet{SVector{N,Complex{BigFloat}}}) where {N} = FE_DirectSolver
default_frame_solver(domain::Domain, basis::FunctionSet{BigFloat}) = FE_DirectSolver
default_frame_solver(domain::Domain, basis::FunctionSet{Complex{BigFloat}}) = FE_DirectSolver
