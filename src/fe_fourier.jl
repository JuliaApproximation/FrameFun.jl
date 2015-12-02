# fe_fourier.jl


function apply!{T,G <: MaskedGrid}(op::Extension, dest, src::DiscreteGridSpace{G}, coef_dest::Array{T}, coef_src::Array{T})
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)

    grid1 = grid(src)
    fill!(coef_dest, zero(T))

    l = 0
    for i in eachindex(grid1.grid)
        if in(i, grid1)
            l = l+1
            coef_dest[i] = coef_src[l]
        end
    end
end


function apply!{T,G <: MaskedGrid}(op::Restriction, dest::DiscreteGridSpace{G}, src, coef_dest::Array{T}, coef_src::Array{T})
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)

    grid1 = grid(dest)

    l = 0
    for i in eachindex(grid1.grid)
        if in(i, grid1)
            l = l+1
            coef_dest[l] = coef_src[i]
        end
    end
end
# This might not be the best way to solve this problem
function discretize_problem{T}(domain::AbstractDomain1d{T}, nt::Tuple{Integer}, tt::Tuple{T}, st::Tuple, basis::DataType, ELT)
    discretize_problem(domain,nt[1],tt[1],st[1],basis,ELT)
end


function discretize_problem{T}(domain::Interval{T}, nt::Int, tt::T, st, basis::DataType, ELT)
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

    rgrid = IndexedSubGrid(grid2, 1, m)

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

    tens_rgrid = TensorProductGrid(ntuple(i->IndexedSubGrid(grid(fbasis2[i]), 1, m[i]), N)...)
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

default_fourier_n(domain::AbstractDomain1d) = 20

default_fourier_n(domain::AbstractDomain2d) = (10, 10)

default_fourier_n(domain::AbstractDomain3d) = (3, 3, 3)

function default_fourier_n{TD,DN,ID,N}(domain::TensorProductDomain{TD,DN,ID,N})
    s=[default_fourier_n(domainlist(domain)[1])...]
    for i=2:ID
        s=[s; default_fourier_n(domainlist(domain)[i])...]
    end
    s=round(Int,s/N)
    tuple(s...)
end

default_fourier_T{T}(domain::AbstractDomain{1,T}) = 2*one(T)
default_fourier_T{N,T}(domain::AbstractDomain{N,T}) = ntuple(i->2*one(T),N)

default_fourier_sampling{T}(domain::AbstractDomain{1,T}) = 2*one(T)
default_fourier_sampling{N,T}(domain::AbstractDomain{N,T}) = ntuple(i->2*one(T),N)


default_fourier_domain_1d() = Interval()

default_fourier_domain_2d() = Circle()

default_fourier_domain_3d() = Sphere()

default_fourier_solver(domain) = FE_ProjectionSolver

#default_fourier_solver(domain::Interval{Float64}) = FE_ProjectionSolver
default_fourier_solver(domain::Interval) = FE_ProjectionSolver

default_fourier_solver(domain::TensorProductDomain) = map(default_fourier_solver,domainlist(domain))    


