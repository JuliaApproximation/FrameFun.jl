# constructors.jl

## """
##   The residual of a SetFun approximation of a Function
## """
## function residual(f::Function, F1::SetFun)
##     F2 = extension_operator(basis(F1))*F1
##     gbasis = BasisFunctions.grid(basis(F2))
##     mask = in(gbasis, domain(F2))
##     Ax = real(full_transform_operator(basis(F2)) * coefficients(F2))[mask]
##     b = sample(gbasis,f)[mask]
##     norm(Ax-b)
## end
## residual(F::SetFun, f::Function) = residual(f,F)

"""
  Create approximation to function with a function set in a domain.

  The number of points is chosen adaptively.
"""
function fun_simple(f::Function, set::FunctionSet, domain::Domain;
    no_checkpoints=200, max_logn_coefs=8, tol=NaN, print_error=false, options...)
  ELT = eltype(f, set)
  N = ndims(set)
  F = nothing
  # TODO Decide which is best
  # tol = default_cutoff(FE_DiscreteProblem(domain, set, 2; options...))
  isequal(tol,NaN) && (tol = 10*10^(4/5*log10(eps(numtype(set)))))
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
  Create approximation to function with with a function set in a domain.

  The number of points is chosen adaptively.
"""
function fun_optimal_N(f::Function, set::FunctionSet, domain::FrameFun.Domain;
    no_checkpoints=200, max_logn_coefs=9, tol=NaN, print_error=false, options...)
  ELT = eltype(f, set)
  N = ndims(set)
  F = nothing
  # TODO Decide which is best
  # tol = default_cutoff(FE_DiscreteProblem(domain, set, 2; options...))
  isequal(tol,NaN) && (tol = 10*10^(4/5*log10(eps(numtype(set)))))
  rgrid = randomgrid(domain, no_checkpoints)
  error = Inf
  random_f=sample(rgrid, f, eltype(f(rgrid[1]...)))
  random_F=zeros(ELT,no_checkpoints)
  N0 = NaN; Nmax = NaN; emax = NaN; Nmin = 1; emin = NaN
  n = 8
  set=resize(set, n)
  while length(set) <= 2^max_logn_coefs
    # println(n, ": ", Nmin, " ", emin, " ", Nmax, " ", emax)
    F=Fun(f, set, domain; options...)
    random_F=F(rgrid)
    error = maximum(abs.(random_F-random_f))
    print_error && (@printf "Error with %d coefficients is %1.3e (%1.3e)\n" (length(set)) error tol)
    if (error<1e-2) && isnan(N0)
            N0 = n
    end
    error < tol ? (Nmax=n; emax=error) : (Nmin=n; emin=error)
    if !isnan(Nmax) && !isnan(Nmin) && (sum(abs.(collect(Nmax)-collect(Nmin))) <= N)
      set=resize(set, Nmax)
      return Fun(f, set, domain; options...)
    end
    # Decide on new n
    if isnan(Nmax)
      n = extension_size(set)
    else
      if isnan(emin)
        return F
      else
        if N==1
          n = round(Int,(log(emin)*Nmin + log(emax)*Nmax)/(log(emin) + log(emax)))
        else
          n = (round.(Int,(log(emin)*collect(Nmin) + log(emax)*collect(Nmax))/(log(emin) + log(emax)))...)
        end
      end
    end
    set=resize(set, n)
  end
  warn("Maximum number of coefficients exceeded, error is $(error)")
  F
end

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
    max_logn_coefs = 7, tol = NaN, options...)
    isequal(tol,NaN) && (tol = 10*10^(4/5*log10(eps(numtype(set)))))
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
struct FunConstructor{N,T}
  space   ::    FunctionSpace{N,T}
  domain  ::    Domain

  FunConstructor{N,T}(space::FunctionSpace, domain::Domain) where {N,T} = new(space, domain)
end

FunConstructor{N,T}(space::FunctionSpace{N,T}, domain::Domain) = FunConstructor{N,T}(space, domain)

domain(constructor::FunConstructor) = constructor.domain
space(constructor::FunConstructor) = constructor.space

# Allows the notation F = FunCSTR(f) with F the approximation and f the function
(constructor::FunConstructor)(f::Function; options...) = construct(constructor, f; options...)

construct(constructor::FunConstructor, f::Function; afun = fun_optimal_N, options...) = afun(f, FunctionSet(space(constructor),0),
    domain(constructor); options...)
