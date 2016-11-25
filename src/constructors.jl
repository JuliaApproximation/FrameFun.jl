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
function fun_simple(f::Function, set::FunctionSet, domain::AbstractDomain;
    no_checkpoints=200, max_logn_coefs=8, tol=NaN, options...)
  ELT = eltype(f, set)
  N = ndims(set)
  F = nothing
  # TODO Decide which is best
  # tol = default_cutoff(FE_DiscreteProblem(domain, set, 2; options...))
  isequal(tol,NaN) && (tol = 10*10^(4/5*log10(eps(numtype(set)))))
  rgrid=random_grid_in_domain(domain,numtype(set);vals=no_checkpoints)
  error = Inf
  random_f=sample(rgrid, f, eltype(f(rgrid[1]...)))
  random_F=zeros(ELT,no_checkpoints)
  n = 8
  set = resize(set,n)
  while length(set) <= 2^(max_logn_coefs)
    set=resize(set,extension_size(set))
    F=Fun(f, set, domain; options...)
    random_F=F(rgrid)
    error = maximum(abs(random_F-random_f))
    @printf "Error with %d coefficients is %1.3e\n" (length(set)) error
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
function fun_optimal_N(f::Function, set::FunctionSet, domain::AbstractDomain;
    no_checkpoints=200, max_logn_coefs=9, tol=NaN, options...)
  ELT = eltype(f, set)
  N = ndims(set)
  F = nothing
  # TODO Decide which is best
  # tol = default_cutoff(FE_DiscreteProblem(domain, set, 2; options...))
  isequal(tol,NaN) && (tol = 10*10^(4/5*log10(eps(numtype(set)))))
  rgrid=random_grid_in_domain(domain,numtype(set);vals=no_checkpoints)
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
    error = maximum(abs(random_F-random_f))
    @printf "Error with %d coefficients is %1.3e (%1.3e)\n" (length(set)) error tol
    if (error<1e-2) && isnan(N0)
      N0 = n
    end
    error < tol ? (Nmax=n; emax=error) : (Nmin=n; emin=error)
    if !isnan(Nmax) && !isnan(Nmin) && sum(abs(collect(Nmax)-collect(Nmin))) <= N
      set=resize(set, Nmax)
      return Fun(f, set, domain; options...)
    end
    # Decide on new n
    if isnan(Nmax)
      n = extension_size(set)
    else
      if N==1
        n = round(Int,(log(emin)*Nmin + log(emax)*Nmax)/(log(emin) + log(emax)))
      else
        n = round(Int,(log(emin)*collect(Nmin) + log(emax)*collect(Nmin))/(log(emin) + log(emax)))
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
  Take a LS of f-p_1, using only phi_2, and call the approximation p_2
  Approximation f-(p_1+p_2) using phi_3,
  ...
"""
function fun_greedy(f::Function, set::FunctionSet1d, domain::AbstractDomain;
    maxn = 200, tol = NaN, options...)
  ELT = eltype(f, set)
  isequal(tol,NaN) && (tol = sqrt(eps(real(ELT))))
  set = promote_eltype(set, ELT)
  coeffs = zeros(ELT,maxn)
  c = zeros(ELT,maxn)
  frame = FrameFuns.domainframe(domain, set)
  frame = resize(frame, maxn)
  g = EquispacedGrid(maxn, left(domain), right(domain))
  fcoeffs = sample(g, f)
  fcoeffs = reshape(fcoeffs,(length(fcoeffs),1))
  n = evaluate_coeffs!(fcoeffs, coeffs, c, set, frame, g, tol; options...)

  F = FrameFun(domain, resize(set, n), coeffs[1:n])
  @printf "Error with %d coefficients is %1.3e (%1.3e)\n" n abserror(f,F) tol
  F
end

function evaluate_coeffs!(fcoeffs::Array, coeffs::Array, c::Array,
    set::FunctionSet1d, frame::DomainFrame, g::AbstractGrid, tol; options...)
  maxn = length(coeffs)
  for n in 1:maxn
    phis = BasisFunctions.evaluation_matrix(resize(frame,n), g)
    c[1:n] = (phis[:,1:n]\fcoeffs)
    coeffs[1:n] += c[1:n]
    fcoeffs -= phis[:,1:n]*c[1:n]
    if minimum(abs(coeffs[1:n]))<tol
      return n
    end

    if n < maxn
        set1 = resize(set,n)
        set2 = resize(set,n+1)
        e = extension_operator(set1,set2)
        c[1:n+1] = apply(extension_operator(set1,set2), coeffs[1:n])
        coeffs[1:n+1] = c[1:n+1]
    end
  end
  return maxn
end

function evaluate_coeffs!(fcoeffs::Array, coeffs::Array, c::Array,
    set::ChebyshevBasis, frame::DomainFrame, g::AbstractGrid, tol; options...)
  maxn = length(coeffs)
  phis = BasisFunctions.evaluation_matrix(frame, g)
  for n in 1:maxn
    c[1:n] = (phis[:,1:n]\fcoeffs)
    coeffs[1:n] += c[1:n]
    fcoeffs -= phis[:,1:n]*c[1:n]
    if norm(coeffs[n])<tol
      return n
    end
  end
  return maxn
end

function evaluate_coeffs!(fcoeffs::Array, coeffs::Array, c::Array,
    set::FourierBasis, frame::DomainFrame, g::AbstractGrid, tol; options...)
  maxn = length(coeffs)
  phis = BasisFunctions.evaluation_matrix(frame, g)
  for n in 1:maxn
    if iseven(n)
      nothing
    else
      I = vcat(maxn-(n>>1)+1:maxn,1:(n>>1)+1)
      c = (phis[:,I]\fcoeffs)
      coeffs[I] += c
      fcoeffs -= phis[:,I]*c
    end
    if norm(coeffs[n])<tol
      return n
    end
  end
  return maxn
end

# Allows following notation
# f2 = cos(f1) with f1 a FrameFun.
for f in (:cos, :sin, :tan, :sinh, :cosh, :tanh,
  :asin, :acos, :atan,
  :sqrt, :cbrt, :exp, :log)
    @eval begin
        Base.$f(F::FrameFun, ; afun = fun_optimal_N, options...) = afun(x->Base.$f(F(x)), basis(F), domain(F); options...)
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
(constructor::FunConstructor)(f::Function; options...) = construct(constructor, f; options...)

construct(constructor::FunConstructor, f::Function; afun = fun_optimal_N, options...) = afun(f, FunctionSet(space(constructor),0),
    domain(constructor); options...)
