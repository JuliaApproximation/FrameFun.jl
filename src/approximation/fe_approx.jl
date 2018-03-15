# fe_fourier.jl


## "An ExpFun is a DictFun based on Fourier series."
## ExpFun(f::Function; n=n, T=T, args...) = Fun(FourierBasis, f ; args...)
## # one of three things should be provided: n(tuple), domain or basis
## # only n

## ExpFun(f::Function, n::Int)
## ExpFun{N}(f::Function, n::Ntuple{N}) = Fun(FourierBasis, f, domain; args...)

## "A ChebyFun is a DictFun based on Chebyshev polynomials."
## ChebyFun(f::Function; args...) = Fun(ChebyshevBasis, f; args...)
## ChebyFun(f::Function, domain; args...) = Fun(ChebyshevBasis, f, domain; args...)


function Fun(f::Function, basis::Dictionary, domain::Domain; options...)
    ELT = codomaintype(f, basis)

    frame = extensionframe(domain, promote_domainsubtype(basis, real(ELT)))
    A = approximation_operator(Span(frame,ELT); options...)
    coef = A * f
    DictFun(domain, dictionary(dest(A)), coef)
end

function fe_problem(basis, domain, sampling_factor; options...)
    frame = ExtensionFrame(domain, basis)
    FE_DiscreteProblem(domain, basis, sampling_factor; options...)
end

function fe_solver(basis, domain; options...)
    frame = ExtensionFrame(domain, basis)
    approximation_operator(Span(frame); options...)
end

# We assume f as a function is type stable.
function codomaintype(f::Function, basis)
    ELT = codomaintype(basis)
    # We only test for the return type in zero
    RT = typeof(f(zero(GeometricSpace{domaintype(basis)})...))
    if (RT <: Complex)
        complex(ELT)
    else
        ELT
    end
end

function oversampled_evaluation_operator(S::Span, D::Domain; sampling_factor=2, incboundary=false, options...)
    B = primaryspan(S)
    # Establish time domain grid
    G, lB = oversampled_grid(D,dictionary(B),sampling_factor)

    op = grid_evaluation_operator(S,gridspace(B,G),G)
    # Add boundary points if necessary
    if incboundary
        BG = boundary(grid(lB), D)
        op = [op; grid_evaluation_operator(S,gridspace(B,BG),BG)]
    end
    (op,scaling_factor(lB))
end

scaling_factor(S::Dictionary) = length(S)
scaling_factor(S::DerivedDict) = scaling_factor(superdict(S))
scaling_factor(S::ChebyshevBasis) = length(S)/2

function discrete_approximation_operator(set::ExtensionSpan; solver = default_frame_solver(domain(set), basisspan(set)), options...)
    (op, scaling) = oversampled_evaluation_operator(basisspan(set),domain(set);options...)
    solver(op; scaling=scaling, options...)
end

primaryspan(span::Span) = span
function primaryspan(span::BasisFunctions.MultiDictSpan)
    elements(span)[findmax(map(length,elements(span)))[2]]
end

struct FE_BestSolver
end

function FE_BestSolver(op::AbstractOperator; scaling=NaN, verbose= false, options...)
    if has_transform(src(op))
        R = estimate_plunge_rank(op)
        if R < size(op, 2)/2
            verbose && println("Estimated plunge rank $R smaller than $(size(op,2))/2 -> projection solver ")
            AZSolver(op; scaling=scaling, verbose=verbose,options...)
        else
            verbose && println("Estimated plunge rank $R greater than $(size(op,2))/2 -> direct solver ")
            DirectSolver(op; verbose=verbose,options...)
        end
    else
        # Don't bother with a fast algorithm if there is no fast transform
        DirectSolver(op, scaling; verbose=verbose, options...)
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

default_frame_solver(domain::Domain, basis::Dictionary{S,SVector{N,BigFloat}}) where {S,N} = DirectSolver
default_frame_solver(domain::Domain, basis::Dictionary{S,SVector{N,Complex{BigFloat}}}) where {S,N} = DirectSolver
default_frame_solver(domain::Domain, basis::Dictionary{BigFloat}) = DirectSolver
default_frame_solver(domain::Domain, basis::Dictionary{Complex{BigFloat}}) = DirectSolver
