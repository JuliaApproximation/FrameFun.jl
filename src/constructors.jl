# constructors.jl

# Create grid consisting of random interior points
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

function AFun(f::Function, set::FunctionSet, domain::AbstractDomain; no_checkpoints=200, max_logn_coefs=8, options...)
  ELT = eltype(f, set)
  # TODO Decide which is best
  # tol = default_cutoff(FE_DiscreteProblem(domain, set, 2; options...))
  tol = 10*10^(4/5*log10(eps(numtype(set))))
  rgrid=random_grid_in_domain(domain,numtype(set);vals=no_checkpoints)

  random_f=sample(rgrid, f, eltype(f(rgrid[1]...)))
  random_F=zeros(ELT,no_checkpoints)
  for logn = 6:max_logn_coefs
    set=resize(set,2^logn)
    F=Fun(f, set, domain; options...)
    random_F=F(rgrid)
    error = reduce(max,abs(random_F-random_f))
    if error<tol
      return F
    end
  end
  warn("Maximum number of coefficients exceeded")
  set=resize(set,2^(max_logn_coefs+1))
  F=Fun(f, set, domain; options...)
end

# function AFun{S<:FunctionSet}(f::Function, ::Type{S}, domain::AbstractDomain,
#     basis_args=(); options...)
#   set = S(16,basis_args...)
#   AFun(f, set, domain; options...)
# end

for f in (:cos, :sin, :tan, :sinh, :cosh, :tanh,
  :asin, :acos, :atan,
  :sqrt, :cbrt, :exp, :log)
    @eval begin
        Base.$f(F::FrameFun) = AFun(x->Base.$f(F(x)),basis(F),domain(F))
    end
end

immutable FunConstructor{N,T}
  space   ::    FunctionSpace{N,T}
  domain  ::    AbstractDomain{N}

  FunConstructor(space::FunctionSpace, domain::AbstractDomain) = new(space, domain)
end

FunConstructor{N,T}(space::FunctionSpace{N,T}, domain::AbstractDomain{N}) = FunConstructor{N,T}(space, domain)

domain(constructor::FunConstructor) = constructor.domain
space(constructor::FunConstructor) = constructor.space
# will not work for multiple set
@compat (constructor::FunConstructor)(f::Function; options...) = AFun(f, FunctionSet(space(constructor),0),
    domain(constructor); options...)

construct(constructor::FunConstructor, f::Function; options...) = AFun(f, FunctionSet(space(constructor),0),
    domain(constructor); options...)
