# constructors.jl
"""
  Create approximation to function with a function set in a domain.

  The number of points is chosen adaptively.
"""
function fun_simple(f::Function, set::FunctionSet, domain::Domain;
    no_checkpoints=200, max_logn_coefs=8, tol=1e-12, print_error=false, options...)
  ELT = rangetype(f, set)
  N = dimension(set)
  F = nothing
  # TODO Decide which is best
  # tol = default_cutoff(FE_DiscreteProblem(domain, set, 2; options...))
  rgrid = randomgrid(domain, no_checkpoints)
  error = Inf
  random_f=sample(rgrid, f, eltype(f(rgrid[1]...)))
  random_F=zeros(ELT,no_checkpoints)
  n = 8
  set = resize(set,n)
  while length(set) <= 2^(max_logn_coefs)
    set=resize(set,extension_size(set))
    F=Fun(f, set, domain; options...)
    random_F=F(rgrid)
    error = maximum(abs.(random_F-random_f))
    print_error && (@printf "Error with %d coefficients is %1.3e\n" (length(set)) error)
    if error<tol
      return F
    end
  end
  warn("Maximum number of coefficients exceeded, error is $(error)")
  F
end
Base.isnan(::Tuple) = false

"""
  Create approximation to function with a function set in a domain.

  The number of points is chosen adaptively.
"""
function fun_optimal_N(f::Function, set::FunctionSet{T}, domain::Domain;
  no_checkpoints=200, max_logn_coefs=8, cutoff=default_cutoff(real(T)), tol=100*cutoff, verbose=false, adaptive_verbose = verbose, return_log=false, randomtest=false, options...) where {T}
  ELT = rangetype(f, set)
  N = dimension(set)
  F = nothing
  rgrid = randomgrid(domain, no_checkpoints)
  error = Inf
  Nmax = NaN; Nmin = 1;
  n = 8
  set=resize(set, n)
  return_log && (log = zeros(T,0,4))
  its = 0
  while length(set) <= 2^max_logn_coefs && its < 100
    # Find new approximation
    F=Fun(f, set, domain; verbose=verbose, cutoff=cutoff, options...)
    # Calculate error
    # Using residual
    error = residual(f, F)/sqrt(length(F))

    # println(n, ": ", Nmin, " ", Nmax, " ", error)
    return_log && (log = [log; n Nmin Nmax error])


    adaptive_verbose  && (@printf "Error with %d coefficients is %1.3e (%1.3e)\n" (length(set)) error tol)

    error < tol ? Nmax=n : Nmin=n

    succeeded_randomtest = true
    # Stop condition
    # Stop if the bounds Nmin and Nmax are determined and if they are close
    if !isnan(Nmax) && !isnan(Nmin) && (sum(abs.(collect(Nmax)-collect(Nmin))) <= N)
        set=resize(set, Nmax)
        F = Fun(f, set, domain; verbose=verbose, cutoff=cutoff, options...)
        #test in 3 random points if randomtest
        randomtestsucceeded = true
        if randomtest
            if maxerror(f, F, vals=3) > tol
                randomtestsucceeded = false
            end
            # return if test succeded
            if randomtestsucceeded
                return_log && (return F, log)
                return F
            else
                # continue if test failed
                adaptive_verbose && (@printf "Random test in 3 points failed\n")
                succeeded_randomtest = false
            end
        else
            return_log && (return F, log)
            return F
        end
    end

    # Decide on new n
    if isnan(Nmax) || !succeeded_randomtest
      # tolerance is not yet met, increase n
      n = extension_size(set)
    else
      # tolerance is reached, find optimal n by subdivision
      if N==1
        # in one dimension Nmin and Nmax are floats
        n = round(Int,(Nmin + Nmax)/2)
      else
        # in multiple dimensions, Nmin and Nmax are tupples
        n = (round.(Int,(collect(Nmin) + collect(Nmax))/2)...)
      end
    end

    # Use new n
    set=resize(set, n)
    its = its + 1
  end
  warn("Maximum number of coefficients exceeded, error is $(error)")
  F
end
default_cutoff(::Type{T}) where T= 10*10^(4/5*log10(eps(real(T))))
default_cutoff(::Type{Float64}) = 1e-16
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
function fun_greedy(f::Function, set::FunctionSet, domain::FrameFun.Domain;
    max_logn_coefs = 7, tol = 1e-12, options...)
    init_n = 4
    set = resize(set,init_n)
    F = Fun(x->0, set, domain; options...)
    for n in init_n:2^max_logn_coefs
        set = resize(set,n)
        p_i = Fun(x->(f(x)-F(x)), set, domain; options...)
        F = F + p_i
        if residual(f, F) < tol
            return F
        end
    end
    F
end

# Allows following notation
# f2 = cos(f1) with f1 a SetFun.
for f in (:cos, :sin, :tan, :sinh, :cosh, :tanh,
  :asin, :acos, :atan,
  :sqrt, :cbrt, :exp, :log)
    @eval begin
        Base.$f(F::SetFun, ; afun = fun_optimal_N, options...) = afun(x->Base.$f(F(x)), basis(F), domain(F); options...)
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

construct(constructor::FunConstructor, f::Function; afun = fun_optimal_N, options...) = afun(f, FunctionSet(space(constructor),0),
    domain(constructor); options...)
