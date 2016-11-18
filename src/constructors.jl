# constructors.jl

"""
  Create grid consisting of random interior points
"""
function random_grid_in_domain{T}(domain::AbstractDomain,::Type{T}=Float64;vals::Int=200)
    # Use the bounding box around the domain
    box = boundingbox(domain)
    N=ndims(domain)
    point=Array{T}(N)
    N == 1 ? points=Array{T}(vals) : points=Array{Vec{N,T}}(vals)
    elements=0
    # Generate some points inside the domain
    while elements < vals
        for j in 1:N
            point[j]=left(box)[j]+(right(box)[j]-left(box)[j])*rand(1)[1]
        end
        N == 1 ? vpoint = point[1] : vpoint = Vec(point...)
        if in(vpoint,domain)
            elements+=1
            points[elements]=vpoint
        end
    end
    return ScatteredGrid(points)
end
"""
  Create approximation to function with with a function set in a domain.

  The number of points is chosen adaptively.
"""
function AFunSimple(f::Function, set::FunctionSet, domain::AbstractDomain;
    no_checkpoints=200, max_logn_coefs=8, tol=NaN, options...)
  ELT = eltype(f, set)
  F = nothing
  # TODO Decide which is best
  # tol = default_cutoff(FE_DiscreteProblem(domain, set, 2; options...))
  isequal(tol,NaN) && (tol = 10*10^(4/5*log10(eps(numtype(set)))))
  rgrid=random_grid_in_domain(domain,numtype(set);vals=no_checkpoints)
  error = Inf
  random_f=sample(rgrid, f, eltype(f(rgrid[1]...)))
  random_F=zeros(ELT,no_checkpoints)
  for logn = 4:max_logn_coefs
    set=resize(set,2^logn)
    F=Fun(f, set, domain; options...)
    random_F=F(rgrid)
    error = maximum(abs(random_F-random_f))
    @printf "Error with %d coefficients is %1.3e\n" (2^(logn*ndims(set))) error
    if error<tol
      return F
    end
  end
  warn("Maximum number of coefficients exceeded, error is $(error)")
  F
end

"""
  Create approximation to function with with a function set in a domain.

  The number of points is chosen adaptively.
"""
function AFun(f::Function, set::FunctionSet, domain::AbstractDomain;
    no_checkpoints=200, max_logn_coefs=9, tol=NaN, options...)
  ELT = eltype(f, set)
  F = nothing
  # TODO Decide which is best
  # tol = default_cutoff(FE_DiscreteProblem(domain, set, 2; options...))
  isequal(tol,NaN) && (tol = 10*10^(4/5*log10(eps(numtype(set)))))
  rgrid=random_grid_in_domain(domain,numtype(set);vals=no_checkpoints)
  error = Inf
  random_f=sample(rgrid, f, eltype(f(rgrid[1]...)))
  random_F=zeros(ELT,no_checkpoints)
  N0 = NaN; Nmax = NaN; emax = NaN; Nmin = 1; emin = NaN
  logn = 4
  while logn <= max_logn_coefs
    if isnan(Nmax)
      n = 2^logn
      logn += 1
    else
      n = round(Int,(log(emin)*Nmin + log(emax)*Nmax)/(log(emin) + log(emax)))
    end
    # println(n, ": ", Nmin, " ", emin, " ", Nmax, " ", emax)
    set=resize(set, n)
    F=Fun(f, set, domain; options...)
    random_F=F(rgrid)
    error = maximum(abs(random_F-random_f))
    @printf "Error with %d coefficients is %1.3e (%1.3e)\n" (n^ndims(F)) error tol
    if (error<1e-2) && isnan(N0)
      N0 = n
    end
    error < tol ? (Nmax=n; emax=error) : (Nmin=n; emin=error)
    if Nmax-Nmin <= 1
      set=resize(set, Nmax)
      return Fun(f, set, domain; options...)
    end
  end
  warn("Maximum number of coefficients exceeded, error is $(error)")
  F
end

"""
  Greedy scheme to get decreasing coefficients.

  Let phi_i i=1..N basisfunctions.
  Take a LS of f with only phi_1, call the approximation p_1
  Take a LS of f-p_1, using only phi_2, and call the approximation p_2
  Approximation f-(p_1+p_2) using phi_3,
  ...
"""
function GreedyFun(f::Function, set::FunctionSet, domain::AbstractDomain;
    maxn = 20, options...)
  ELT = eltype(f, set)
  coeffs = zeros(ELT,maxn)
  GreedyFun(f, promote_eltype(set, ELT), domain, coeffs; options...)
end

# TODO assumes grid of set exists
function GreedyFun(f::Function, set::FunctionSet, domain::AbstractDomain, coeffs;
    tol=NaN, options...)
  maxn = length(coeffs)
  frame = FrameFuns.domainframe(domain, set)
  frame = resize(frame, maxn)
  fcoeffs = sample(grid(frame), f)
  println(ctranspose(fcoeffs))
  for n in 1:maxn
    phi = BasisFunctions.evaluation_matrix(frame[n], BasisFunctions.grid(frame))[:]
    coeffs[n] = (phi\fcoeffs)[1]
    p = coeffs[n]*phi
    fcoeffs -= p
  end
  FrameFun(domain, resize(set, length(coeffs)), coeffs)
end


# Allows following notation
# f2 = cos(f1) with f1 a FrameFun.
for f in (:cos, :sin, :tan, :sinh, :cosh, :tanh,
  :asin, :acos, :atan,
  :sqrt, :cbrt, :exp, :log)
    @eval begin
        Base.$f(F::FrameFun; options...) = AFun(x->Base.$f(F(x)), basis(F), domain(F); options...)
    end
end

# A FunConstructor approximates a function in a domain given a (function)space
immutable FunConstructor{N,T}
  space   ::    FunctionSpace{N,T}
  domain  ::    AbstractDomain{N}

  FunConstructor(space::FunctionSpace, domain::AbstractDomain) = new(space, domain)
end

FunConstructor{N,T}(space::FunctionSpace{N,T}, domain::AbstractDomain{N}) = FunConstructor{N,T}(space, domain)

domain(constructor::FunConstructor) = constructor.domain
space(constructor::FunConstructor) = constructor.space

# Allows the notation F = FunCSTR(f) with F the approximation and f the function
(constructor::FunConstructor)(f::Function; options...) = AFun(f, FunctionSet(space(constructor),0),
    domain(constructor); options...)

construct(constructor::FunConstructor, f::Function; options...) = AFun(f, FunctionSet(space(constructor),0),
    domain(constructor); options...)
