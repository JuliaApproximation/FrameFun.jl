# fe_fourier.jl


# This might not be the best way to solve this problem
function discretize_problem{T}(domain::AbstractDomain1d{T}, nt::Tuple{Integer}, tt::Tuple{T}, st::Tuple, basis::DataType, ELT)
    discretize_problem(domain,nt[1],tt[1],st[1],basis,ELT)
end


function discretize_problem{T}(domain::Interval{T}, nt::Int, tt, st, basis::Type{FourierBasis}, ELT)
    n = 2*nt+1
    m = 2*round(Int, nt.*st)+1
    t = (tt.*(m-1)/2).*(2./(m-1))
    l = round(Int, t.*(m-1))

    t = (l*one(T)) / ((m-1)*one(T))
    a = left(domain)
    b = right(domain)

    fbasis1 = basis(n, a, b + (b-a)*(t-1))
    fbasis2 = basis(l, a, b + (b-a)*(t-1))

    grid1 = grid(fbasis1)
    grid2 = grid(fbasis2)

    rgrid = IndexSubGrid(grid2, 1, m)
    
    tbasis1 = DiscreteGridSpace(grid1, ELT)
    tbasis2 = DiscreteGridSpace(grid2, ELT)

    tbasis_restricted = DiscreteGridSpace(rgrid, ELT)

    FE_DiscreteProblem(domain, fbasis1, fbasis2, tbasis1, tbasis2, tbasis_restricted)
end


function discretize_problem{T}(domain::AbstractDomain1d{T}, nt::Int, tt, st, basis::DataType, ELT)
    n = 2*nt+1
    m = 2*round(Int, nt.*st)+1
    t = convert(numtype(domain),(tt.*(m-1)/2).*(2./(m-1)))
    l = round(Int, t.*(m-1))

    t = (l*one(T)) / ((m-1)*one(T))

    a = left(domain)
    b = right(domain)

    fbasis1 = basis(n, a, b + (b-a)*(t-1))
    fbasis2 = basis(l, a, b + (b-a)*(t-1))

    grid1 = grid(fbasis1)
    grid2 = grid(fbasis2)

    rgrid = MaskedGrid(grid2, domain)
    
    tbasis1 = DiscreteGridSpace(grid1, ELT)
    tbasis2 = DiscreteGridSpace(grid2, ELT)
    
    tbasis_restricted = DiscreteGridSpace(rgrid, ELT)

    FE_DiscreteProblem(domain, fbasis1, fbasis2, tbasis1, tbasis2, tbasis_restricted)
    
end

function discretize_problem{N,T}(domain::AbstractDomain{N,T}, nt::Tuple, tt::Tuple, st::Tuple, basis::DataType, ELT)
    n = 2*[nt...]+1
    m = 2*round(Int, [nt...].*[st...])+1
    tt = round(Int,[tt...].*(m-1)/2).*(2./(m-1))
    l = round(Int, tt.*(m-1))
    fbasis1=Array{basis}(N)
    fbasis2=Array{basis}(N)
    bbox = box(domain)
    for i=1:N
        t = (l[1]*one(T))/((m[1]-1)*one(T))
        fbasis1[i] = basis(n[i], left(bbox)[i], right(bbox)[i] + (right(bbox)[i]-left(bbox)[i])*(t-1))
        fbasis2[i] = basis(l[i], left(bbox)[i], right(bbox)[i] + (right(bbox)[i]-left(bbox)[i])*(t-1))
    end
    tens_fbasis1 = TensorProductSet(fbasis1...)
    tens_fbasis2 = TensorProductSet(fbasis2...)

    tens_grid1 = TensorProductGrid(map(x->grid(x),fbasis1)...)
    tens_grid2 = TensorProductGrid(map(x->grid(x),fbasis2)...)

    tens_rgrid = TensorProductGrid(ntuple(i->IndexSubGrid(grid(fbasis2[i]), 1, m[i]), N)...)
    tens_tbasis1 = TensorProductSet(map(x->DiscreteGridSpace(grid(x),ELT), fbasis1)...)
    tens_tbasis2 = TensorProductSet(map(x->DiscreteGridSpace(grid(x),ELT), fbasis2)...)

    tbasis_restricted = DiscreteGridSpace(MaskedGrid(tens_grid2, domain),ELT)

    FE_DiscreteProblem(domain, tens_fbasis1, tens_fbasis2, tens_tbasis1, tens_tbasis2, tbasis_restricted)
end

function discretize_problem{TD,DN,ID,N,T}(domain::TensorProductDomain{TD,DN,ID,N,T}, nt::Tuple, tt::Tuple, st::Tuple, basis::DataType, ELT)
    problems=FE_Problem[]
    dc = 1
    for i=1:ID
        push!(problems,discretize_problem(subdomain(domain,i),nt[dc:dc+DN[i]-1],tt[dc:dc+DN[i]-1],st[dc:dc+DN[i]-1],basis,ELT))
        dc=dc+DN[i]
    end
    problem = FE_TensorProductProblem(problems...)
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


# The default basis is a FourierBasis
default_frame_domain_1d() = default_frame_domain_1d(FourierBasis)
default_frame_domain_2d() = default_frame_domain_2d(FourierBasis)
default_frame_domain_3d() = default_frame_domain_3d(FourierBasis)
default_frame_n(domain) = default_frame_n(domain, FourierBasis)
default_frame_T(domain) = default_frame_T(domain, FourierBasis)
default_frame_sampling(domain) = default_frame_sampling(domain, FourierBasis)
default_frame_solver(domain) = default_frame_solver(domain, FourierBasis)

default_frame_domain_1d{Basis}(::Type{Basis}) = Interval()
default_frame_domain_2d{Basis}(::Type{Basis}) = Circle()
default_frame_domain_3d{Basis}(::Type{Basis}) = Sphere()


default_frame_n{Basis}(domain::AbstractDomain1d, ::Type{Basis}) = 20
default_frame_n{Basis}(domain::AbstractDomain2d, ::Type{Basis}) = (10, 10)
default_frame_n{Basis}(domain::AbstractDomain3d, ::Type{Basis}) = (3, 3, 3)

default_frame_n(domain::AbstractDomain1d, ::Type{FourierBasis}) = 21


function default_frame_n{Basis,TD,DN,ID,N}(domain::TensorProductDomain{TD,DN,ID,N}, ::Type{Basis})
    s = [default_frame_n(domainlist(domain)[1], Basis)...]
    for i = 2:ID
        s=[s; default_frame_n(domainlist(domain)[i], Basis)...]
    end
    s = round(Int,s/N)
    tuple(s...)
end

default_frame_T{T,Basis <: FunctionSet}(domain::AbstractDomain{1,T}, ::Type{Basis}) = 2*one(T)
default_frame_T{N,T,Basis <: FunctionSet}(domain::AbstractDomain{N,T}, ::Type{Basis}) = ntuple(i->2*one(T),N)


default_frame_sampling{T}(domain::AbstractDomain{1,T}, Basis) = 2*one(T)
default_frame_sampling{N,T}(domain::AbstractDomain{N,T}, Basis) = ntuple(i->2*one(T),N)


default_frame_solver{Basis}(domain, ::Type{Basis}) = FE_ProjectionSolver

#default_frame_solver(domain::Interval{Float64}) = FE_ProjectionSolver
default_frame_solver{Basis}(domain::Interval, ::Type{Basis}) = FE_ProjectionSolver

default_frame_solver{Basis}(domain::TensorProductDomain, ::Type{Basis}) = map(default_frame_solver,domainlist(domain))


