# fe_fourier.jl


"An ExpFun is a FrameFun based on Fourier series."
ExpFun(f::Function; args...) = Fun(FourierBasis, f; args...)
ExpFun(f::Function, domain; args...) = Fun(FourierBasis, f, domain; args...)

"A ChebyFun is a FrameFun based on Chebyshev polynomials."
ChebyFun(f::Function; args...) = Fun(ChebyshevBasis, f; args...)
ChebyFun(f::Function, domain; args...) = Fun(ChebyshevBasis, f, domain; args...)


function Fun{Basis <: FunctionSet}(::Type{Basis}, f::Function, domain = default_frame_domain_1d(); args...)
    ELT = eltype(f, domain, Basis)
    (problem,solver) = fe_problem(domain, Basis, ELT; args...)
    coef = solve(solver, f, problem)
    FrameFun(domain, frequency_basis(problem), coef)
end


function eltype{Basis <: AbstractBasis}(f::Function, domain, ::Type{Basis})
    ELT = numtype(domain)
    RT = Base.return_types(f,fill(numtype(domain),dim(domain)))
    if length(RT) > 0
        if isreal(Basis) == False || (RT[1] <: Complex)
            Complex{ELT}
        else
            ELT
        end
    else
        if isreal(Basis) == False
            Complex{ELT}
        else
            ELT
        end
    end
end



"""
Construct an FE problem for the given domain, using default values if necessary.
"""
function fe_problem(domain, Basis, ELT;
    n = default_frame_n(domain, Basis),
    T = default_frame_T(domain, Basis),
    s = default_frame_sampling(domain, Basis),
    solver = default_frame_solver(domain),
    args...)
    
    problem = discretize_problem(domain, n, T, s, Basis, ELT)
    sol = solver(problem)

    (problem, sol)
end





# This might not be the best way to solve this problem
function discretize_problem(domain::AbstractDomain1d, nt::Tuple, tt::Tuple, st::Tuple, Basis, ELT)
    discretize_problem(domain, nt[1], tt[1], st[1], Basis, ELT)
end


# This routine allows to generalize discretize_problem, but it is very specific.
# TODO: make this more general. The parameter 'm' is rather ugly.
suitable_subgrid(grid, domain, basis, m) = MaskedGrid(grid, domain)

suitable_subgrid(grid, domain::Interval, ::Type{FourierBasis}, m) = IndexSubGrid(grid, 1, m)


# We should be doing this adaptively, in order to guarantee a certain oversampling rate
function discretize_problem{T}(domain::AbstractDomain1d{T}, n::Int, tt, st, Basis, ELT)
    m = round(Int, n.*st)+1
    t = convert(numtype(domain),(tt.*(m-1)/2).*(2./(m-1)))
    l = round(Int, t.*(m-1))

    t = (l*one(T)) / ((m-1)*one(T))

    a = left(domain)
    b = right(domain)

    fbasis1 = Basis(n, a, b + (b-a)*(t-1))
    fbasis2 = Basis(l, a, b + (b-a)*(t-1))

    grid1 = grid(fbasis1)
    grid2 = grid(fbasis2)

    # For FourierBasis and Interval:
    #rgrid = IndexSubGrid(grid2, 1, m)
    # For anything else:
    #rgrid = MaskedGrid(grid2, domain)
    rgrid = suitable_subgrid(grid2, domain, Basis, m)
    
    tbasis1 = DiscreteGridSpace(grid1, ELT)
    tbasis2 = DiscreteGridSpace(grid2, ELT)
    
    tbasis_restricted = DiscreteGridSpace(rgrid, ELT)

    FE_DiscreteProblem(domain, fbasis1, fbasis2, tbasis1, tbasis2, tbasis_restricted)
    
end

function discretize_problem{N,T}(domain::AbstractDomain{N,T}, nt::Tuple, tt::Tuple, st::Tuple, Basis, ELT)
    n = [nt...]
    m = round(Int, [nt...].*[st...])+1
    tt = round(Int,[tt...].*(m-1)/2).*(2./(m-1))
    l = round(Int, tt.*(m-1))
    fbasis1 = Array{Basis}(N)
    fbasis2 = Array{Basis}(N)
    bbox = box(domain)
    for i=1:N
        t = (l[1]*one(T))/((m[1]-1)*one(T))
        fbasis1[i] = Basis(n[i], left(bbox)[i], right(bbox)[i] + (right(bbox)[i]-left(bbox)[i])*(t-1))
        fbasis2[i] = Basis(l[i], left(bbox)[i], right(bbox)[i] + (right(bbox)[i]-left(bbox)[i])*(t-1))
    end
    tens_fbasis1 = TensorProductSet(fbasis1...)
    tens_fbasis2 = TensorProductSet(fbasis2...)

    tens_grid1 = TensorProductGrid(map(grid, fbasis1)...)
    tens_grid2 = TensorProductGrid(map(grid, fbasis2)...)

    tens_rgrid = TensorProductGrid(ntuple(i->IndexSubGrid(grid(fbasis2[i]), 1, m[i]), N)...)
    tens_tbasis1 = TensorProductSet(map(x->DiscreteGridSpace(grid(x),ELT), fbasis1)...)
    tens_tbasis2 = TensorProductSet(map(x->DiscreteGridSpace(grid(x),ELT), fbasis2)...)

    tbasis_restricted = DiscreteGridSpace(MaskedGrid(tens_grid2, domain), ELT)

    FE_DiscreteProblem(domain, tens_fbasis1, tens_fbasis2, tens_tbasis1, tens_tbasis2, tbasis_restricted)
end

function discretize_problem{TD,DN,LEN,N,T}(domain::TensorProductDomain{TD,DN,LEN,N,T}, nt::Tuple, tt::Tuple, st::Tuple, Basis, ELT)
    problems = FE_Problem[]
    dc = 1
    for i = 1:LEN
        range = dc:dc+DN[i]-1
        p = discretize_problem(subdomain(domain,i), nt[range], tt[range], st[range], Basis, ELT)
        push!(problems, p)
        dc += DN[i]
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
for op in (:default_frame_1d, :default_frame_2d, :default_frame_3d)
    @eval $op() = $op(FourierBasis)
end

for op in (:default_frame_n, :default_frame_T, :default_frame_sampling, :default_frame_solver)
    @eval $op(domain) = $op(domain, FourierBasis)
end

default_frame_domain_1d(Basis) = Interval()
default_frame_domain_2d(Basis) = Disk()
default_frame_domain_3d(Basis) = Ball()


default_frame_n(domain::AbstractDomain1d, Basis) = 61
default_frame_n(domain::AbstractDomain2d, Basis) = (31, 31)
default_frame_n(domain::AbstractDomain3d, Basis) = (7, 7, 7)


function default_frame_n(domain::TensorProductDomain, Basis)
    s = [default_frame_n(domainlist(domain)[1], Basis)...]
    for i = 2:tp_length(domain)
        s = [s; default_frame_n(domainlist(domain)[i], Basis)...]
    end
    s = round(Int,s/dim(domain))
    tuple(s...)
end

default_frame_T{T}(domain::AbstractDomain{1,T}, Basis) = 2*one(T)
default_frame_T{N,T}(domain::AbstractDomain{N,T}, Basis) = ntuple(i->2*one(T),N)


default_frame_sampling{T}(domain::AbstractDomain{1,T}, Basis) = 2*one(T)
default_frame_sampling{N,T}(domain::AbstractDomain{N,T}, Basis) = ntuple(i->2*one(T),N)


default_frame_solver(domain, Basis) = FE_ProjectionSolver

default_frame_solver{N}(domain::AbstractDomain{N,BigFloat}, Basis) = FE_DirectSolver

# Does not seem to work (anymore)
#default_frame_solver(domain::TensorProductDomain, Basis) = map(default_frame_solver, domainlist(domain))

