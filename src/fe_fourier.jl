# fe_fourier.jl



# Perform the copy without additional allocation
@generated function copy_ranges!{N}(c, d, c_ranges::NTuple{N}, d_ranges::NTuple{N})
    quote
        @nloops $N i x->1:length(c_ranges[x]) begin
         (@nref $N d x->d_ranges[x][i_x]) = (@nref $N c x->c_ranges[x][i_x])
        end
    end
end


## This code is in anticipation of a working TensorProductOperator implementation.
# TODO: fix Chebyshev transform and remove these functions below

apply_tensor!(op::Extension, dest::TensorProductSet, src::TensorProductSet, dest1::ChebyshevBasis, src1::ChebyshevBasis, coef_dest, coef_src) = 
    reshape_N_to_L_cheby!(coef_src, coef_dest, size(coef_src), size(coef_dest))

apply_tensor!(op::Restriction, dest::TensorProductSet, src::TensorProductSet, dest1::ChebyshevBasis, src1::ChebyshevBasis, coef_dest, coef_src) =
    reshape_L_to_N_cheby!(coef_dest, coef_src, size(coef_dest), size(coef_src))

# Reshape functions: we want to efficiently copy the data from a vector of length N to a larger vector of length L.
# Hard to do with Cartesian, but it can be done for any dimension recursively.
# TODO: find a better way!
function reshape_N_to_L_cheby!{N}(c, d, n::NTuple{N}, l::NTuple{N})
    nh = map(x->div(x-1,2), n)
    # First zero out d
    fill!(d, 0)
    reshape_N_to_L_rec_cheby!(c, d, (), (), nh, n, l)
end

function reshape_N_to_L_rec_cheby!{N}(c, d, c_ranges, d_ranges, nh::NTuple{N}, n, l)
    reshape_N_to_L_rec_cheby!(c, d, tuple(c_ranges...,1:n[1]), tuple(d_ranges...,1:n[1]), nh[2:end], n[2:end], l[2:end])
end

# The end of the recursion: perform the actual copy
function reshape_N_to_L_rec_cheby!(c, d, c_ranges, d_ranges, nh::NTuple{1}, n::NTuple{1}, l::NTuple{1})
    # Currently, the two lines below do some allocation. Using views is not a great improvement.
    # d[d_ranges...,1:nh[1]+1] = c[c_ranges...,1:nh[1]+1]
    # d[d_ranges...,l[1]-nh[1]+1:l[1]] = c[c_ranges...,n[1]-nh[1]+1:n[1]]
    copy_ranges!(c, d, tuple(c_ranges...,1:n[1]), tuple(d_ranges...,1:n[1]))
end


function reshape_L_to_N_cheby!{N}(c, d, n::NTuple{N}, l::NTuple{N})
    nh = map(x->div(x-1, 2), n)
    reshape_L_to_N_rec_cheby!(c, d, (), (), nh, n, l)
end

function reshape_L_to_N_rec_cheby!{N}(c, d, c_ranges, d_ranges, nh::NTuple{N}, n, l)
    reshape_L_to_N_rec_cheby!(c, d, tuple(c_ranges...,1:n[1]), tuple(d_ranges...,1:n[1]), nh[2:end], n[2:end], l[2:end])
end

# The end of the recursion: perform the actual copy
function reshape_L_to_N_rec_cheby!(c, d, c_ranges, d_ranges, nh::NTuple{1}, n::NTuple{1}, l::NTuple{1})
    copy_ranges!(d, c, tuple(d_ranges...,1:n[1]), tuple(c_ranges...,1:n[1]))
end






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


function discretize_problem{T}(domain::AbstractDomain1d{T}, nt::Int, tt::T, st, basis::DataType, ELT)
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

function discretize_problem{T}(domain::AbstractDomain1d{T}, n::Tuple, t::Tuple, s::Tuple)
    discretize_problem(domain,n[1],t[1],s[1])
end

default_fourier_domain_1d() = Interval()

default_fourier_domain_2d() = Circle()

default_fourier_domain_3d() = Sphere()

default_fourier_solver(domain) = FE_ProjectionSolver

#default_fourier_solver(domain::Interval{Float64}) = FE_ProjectionSolver
default_fourier_solver(domain::Interval) = FE_ProjectionSolver

default_fourier_solver(domain::TensorProductDomain) = map(default_fourier_solver,domainlist(domain))    


