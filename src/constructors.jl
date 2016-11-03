# constructors.jl

# Create grid consisting of random interior points
function random_grid_in_domain{T}(domain::AbstractDomain,::Type{T}=Float64;vals::Int=200)
    # Use the bounding box around the domain
    box = boundingbox(domain)
    N=ndims(domain)
    point=Array{T}(N)
    N == 1 ? points=Array{T}(vals) : points=Array{Vec{N,T}}(vals)
    elements=0
    # Generate some points inside the domain, and compare with the target function
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

function Fun{S<:FunctionSet}(f::Function, ::Type{S}, domain::AbstractDomain, basis_args=(); no_checkpoints=200,options...)
  set=S(1,basis_args...)
  ELT = eltype(f, set)
  tol = 10*10^(4/5*log10(eps(numtype(set))))
  rgrid=random_grid_in_domain(domain,numtype(set);vals=no_checkpoints)
  random_f=sample(rgrid, f)
  random_F=zeros(ELT,no_checkpoints)
  for logn = 6:8
    set=S(2^logn,basis_args...)
    F=Fun(f, set, domain; options...)
    random_F=F(rgrid)
    error = reduce(max,abs(random_F-random_f))
    if error<tol
      return F
    end
  end
  warn("Maximum number of coefficients reached exceeded")
  set=S(2^9,basis_args...)
  F=Fun(f, set, domain; options...)
end
