
function _trapnorm(f::Function, D::ExtensionFrame, oversamplingfactor)
    g = FrameFun.oversampled_grid(domain(D), basis(D), oversamplingfactor = oversamplingfactor)[1]
    sqrt(stepsize(g)*norm([f(x) for x in g])^2)
end

BasisFunctions.stepsize(s::IndexSubGrid) = BasisFunctions.stepsize(supergrid(s))

random_test(f::Function, F::DictFun, tolerance, no_samples::Int) =
    abserror(f,F;vals=no_samples) < tolerance

abs_coefficient_test(f::Function, F::DictFun, abscoef, oversamplingfactor) =
    norm(coefficients(F)) < abscoef*_trapnorm(f, dictionary(F), oversamplingfactor)

rel_coefficient_test(f::Function, F::DictFun, relcoef, prev_cnorm, tolerance) =
    abs(norm(coefficients(F))-prev_cnorm) < relcoef*tolerance

function coefficient_test(f::Function, F::DictFun; abscoef=nothing, relcoef=nothing, tolerance=nothing, oversamplingfactor=1, prev_cnorm=nothing, options...)
    if relcoef != nothing
        rel_coefficient_test(f, F, relcoef, prev_cnorm, tolerance)
    elseif abscoef != nothing
        abs_coefficient_test(f, F, abscoef, oversamplingfactor)
    else
        error()
    end
end



"""
  Create approximation to function with a given dictionary in a domain.

  The number of points is chosen adaptively, and guarantees the tolerance, but may not be optimal.
"""
function fun_simple(f::Function, dict::Dictionary, domain::Domain;
        no_checkpoints=3, max_logn_coefs=8, tol=1e-12, abscoef=nothing, verbose=false, adaptive_verbose = verbose, options...)
    ELT = codomaintype(f, dict)
    N = dimension(dict)
    F = nothing
    rgrid = randomgrid(domain, no_checkpoints)
    error = Inf
    random_f=sample(rgrid, f, eltype(f(rgrid[1]...)))
    random_F=zeros(ELT,no_checkpoints)
    n = 8
    dict = resize(dict,n)
    while length(dict) <= 1<<(max_logn_coefs)
        dict = resize(dict,extension_size(dict))
        F = Fun(f, dict, domain; verbose=verbose, options...)
        random_F = F(rgrid)
        error = maximum(abs.(random_F-random_f))
        adaptive_verbose && (@printf "Error with %d coefficients is %1.3e\n" (length(dict)) error)
        (adaptive_verbose & (abscoef != nothing))  && (@printf "Norm with %d coefficients is %1.3e (%1.3e)\n" (length(dict)) norm(coefficients(F)) abscoef*_trapnorm(f, dictionary(F), 1))
        if (error < tol) && ((abscoef == nothing) || abs_coefficient_test(f, F, abscoef, 1))
            return F
        end
    end
    @warn("Maximum number of coefficients exceeded, error is $(error)")
    F
end

Base.isnan(::Tuple) = false

are_close(N1, N2) =
    !isnan(N1) && !isnan(N2) && (sum(abs.(collect(N1)-collect(N2))) <= length(N1))

"""
  Create approximation to function with a dictionary in a domain.

  The number of points is chosen adaptively and optimally.
"""
function fun_optimal_N(f::Function, dict::Dictionary{S,T}, domain::Domain;
        no_checkpoints=3, max_logn_coefs=8, threshold=default_threshold(real(S)), tol=100*threshold, verbose=false, adaptive_verbose = verbose, return_log=false, randomtest=false, abscoef=nothing, relcoef=nothing, options...) where {S,T}
    coefficienttest = (nothing!=abscoef) | (nothing!=relcoef)

    ELT = codomaintype(f, dict)
    N = dimension(dict)
    F = nothing
    rgrid = randomgrid(domain, no_checkpoints)
    Nmax = NaN; Nmin = 1;
    n = 8
    dict=resize(dict, n)
    return_log && (log = zeros(T,0,4))
    its = 0
    error=-1
    while length(dict) <= 2^max_logn_coefs && its < 100
        # Find new approximation
        F=Fun(f, dict, domain; verbose=verbose, threshold=threshold, options...)

        # Using residual
        error = residual(f, F)/sqrt(length(F))
        # println(n, ": ", Nmin, " ", Nmax, " ", error)
        return_log && (log = [log; n Nmin Nmax error])
        adaptive_verbose  && (@printf "Error with %d coefficients is %1.3e (%1.3e)\n" (length(dict)) error tol)

        errorsucces = (error < tol)

        randomsuccess = randomtest ?
            random_test(f, F, tol, 3) : true
        (adaptive_verbose&randomtest)  && (randomsuccess ?
            (@printf "Random test in 3 points succeeded\n") : (@printf "Random test in 3 points failed\n"))

        coefsuccess = coefficienttest ?
            coefficient_test(f, F; tolerance=tol, abscoef=abscoef, relcoef=relcoef, options...) : true
        (adaptive_verbose&coefficienttest)  && (coefsuccess ?
            (@printf "Coefficient test succeeded %1.3e\n" norm(coefficients(F))) :
            (@printf "Coefficient test failed %1.3e\n" norm(coefficients(F))))

        success = broadcast(&, errorsucces, coefsuccess, randomsuccess)
        success ? Nmax=n : Nmin=n

        # Stop condition
        # If the bounds Nmin and Nmax are determined and if they are close
        if are_close(Nmin, Nmax)
            dict=resize(dict, Nmax)
            F = Fun(f, dict, domain; verbose=verbose, threshold=threshold, options...)

            return_log && (return F, log)
            return F
        end

        # Decide on new n
        if isnan(Nmax)
            # tolerance is not yet met, increase n
            n = extension_size(dict)
        else
            # tolerance is reached, find optimal n by subdivision
            if N==1
                # in one dimension Nmin and Nmax are floats
                n = round(Int,(Nmin + Nmax)/2)
            else
                # in multiple dimensions, Nmin and Nmax are tupples
                n = (round.(Int,(collect(Nmin) + collect(Nmax))/2)...,)
            end
        end

        # Use new n
        dict=resize(dict, n)
        its = its + 1
    end
    @warn("Maximum number of coefficients exceeded, error is $(error)")
    F
end


default_threshold(::Type{T}) where T= 10*10^(4/5*log10(eps(real(T))))
default_threshold(::Type{Float64}) = 1e-16
Base.real(::Type{Tuple{N,T}}) where {N,T} = real(T)

"""
  Greedy scheme to get decreasing coefficients.

  Let phi_i i=1..N basisfunctions.
  Take a LS of f with only phi_1, call the approximation p_1
  Take a LS of f-p_1, using only {phi_1,phi_2}, and call the approximation p_2
  Approximation f-(p_1+p_2) using {phi_k}_{k=1}^3,
  ...
  Stop the iteration when the residu of the approximation is smaller than the tolerance
  or at the maximum number of iterations.
"""
function fun_greedy(f::Function, dict::Dictionary, domain::FrameFun.Domain;
    max_logn_coefs = 7, tol = 1e-12, options...)
    init_n = 4
    dict = resize(dict,init_n)
    F = Fun(x->0, dict, domain; options...)
    for n in init_n:2^max_logn_coefs
        dict = resize(dict,n)
        p_i = Fun(x->(f(x)-F(x)), dict, domain; options...)
        F = F + p_i
        if residual(f, F) < tol
            return F
        end
    end
    F
end

function FourierFun(f::Function, ELT = Float64; T=ELT(2), Omega=Interval(-ELT(-1),ELT(1)), Gamma=T*Omega, options...)
    B = FourierBasis(1, leftendpoint(Omega), rightendpoint(Omega), ELT)
    D = Omega
    frame = extensionframe(B, Gamma)
end

# Allows following notation
# f2 = cos(f1) with f1 a DictFun.
for f in (:cos, :sin, :tan, :sinh, :cosh, :tanh,
    :asin, :acos, :atan,
    :sqrt, :cbrt, :exp, :log)
    @eval begin
        Base.$f(F::DictFun, ; afun = fun_optimal_N, options...) = afun(x->Base.$f(F(x)), basis(F), domain(F); options...)
    end
end

# A FunConstructor approximates a function in a domain given a (function)space
struct FunConstructor{T}
    space   ::    FunctionSpace{T}
    domain  ::    Domain

    FunConstructor{T}(space::FunctionSpace, domain::Domain) where {T} = new(space, domain)
end

FunConstructor(space::FunctionSpace{T}, domain::Domain) where {T} = FunConstructor{T}(space, domain)

domain(constructor::FunConstructor) = constructor.domain
space(constructor::FunConstructor) = constructor.space

# Allows the notation F = FunCSTR(f) with F the approximation and f the function
(constructor::FunConstructor)(f::Function; options...) = construct(constructor, f; options...)

construct(constructor::FunConstructor, f::Function; afun = fun_optimal_N, options...) = afun(f, Dictionary(space(constructor),0),
    domain(constructor); options...)
