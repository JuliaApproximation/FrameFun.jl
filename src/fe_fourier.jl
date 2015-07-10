# fe_fourier.jl



# Reshape functions: we want to efficiently copy the data from a vector of length N to a larger vector of length L.
# Hard to do with Cartesian, but it can be done for any dimension recursively.
# TODO: find a better way!
function reshape_N_to_L!{N}(c, d, n::NTuple{N}, l::NTuple{N})
    nh = map(x->div(x-1,2), n)
    # First zero out d
    fill!(d, 0)
    reshape_N_to_L_rec!(c, d, (), (), nh, n, l)
    reshape_N_to_L_rec!(c, d, (), (), nh, n, l)
end

function reshape_N_to_L_rec!{N}(c, d, c_ranges, d_ranges, nh::NTuple{N}, n, l)
    reshape_N_to_L_rec!(c, d, tuple(c_ranges...,1:nh[1]+1), tuple(d_ranges...,1:nh[1]+1), nh[2:end], n[2:end], l[2:end])
    reshape_N_to_L_rec!(c, d, tuple(c_ranges...,n[1]-nh[1]+1:n[1]), tuple(d_ranges...,l[1]-nh[1]+1:l[1]), nh[2:end], n[2:end], l[2:end])
end

# The end of the recursion: perform the actual copy
function reshape_N_to_L_rec!(c, d, c_ranges, d_ranges, nh::NTuple{1}, n::NTuple{1}, l::NTuple{1})
    # Currently, the two lines below do some allocation. Using views is not a great improvement.
    # d[d_ranges...,1:nh[1]+1] = c[c_ranges...,1:nh[1]+1]
    # d[d_ranges...,l[1]-nh[1]+1:l[1]] = c[c_ranges...,n[1]-nh[1]+1:n[1]]
    copy_ranges!(c, d, tuple(c_ranges...,1:nh[1]+1), tuple(d_ranges...,1:nh[1]+1))
    copy_ranges!(c, d, tuple(c_ranges...,n[1]-nh[1]+1:n[1]), tuple(d_ranges...,l[1]-nh[1]+1:l[1]))
end


function reshape_L_to_N!{N}(c, d, n::NTuple{N}, l::NTuple{N})
    nh = map(x->div(x-1, 2), n)
    reshape_L_to_N_rec!(c, d, (), (), nh, n, l)
    reshape_L_to_N_rec!(c, d, (), (), nh, n, l)
end

function reshape_L_to_N_rec!{N}(c, d, c_ranges, d_ranges, nh::NTuple{N}, n, l)
    reshape_L_to_N_rec!(c, d, tuple(c_ranges...,1:nh[1]+1), tuple(d_ranges...,1:nh[1]+1), nh[2:end], n[2:end], l[2:end])
    reshape_L_to_N_rec!(c, d, tuple(c_ranges...,n[1]-nh[1]+1:n[1]), tuple(d_ranges...,l[1]-nh[1]+1:l[1]), nh[2:end], n[2:end], l[2:end])
end

# The end of the recursion: perform the actual copy
function reshape_L_to_N_rec!(c, d, c_ranges, d_ranges, nh::NTuple{1}, n::NTuple{1}, l::NTuple{1})
    copy_ranges!(d, c, tuple(d_ranges...,1:nh[1]+1), tuple(c_ranges...,1:nh[1]+1))
    copy_ranges!(d, c, tuple(d_ranges...,l[1]-nh[1]+1:l[1]), tuple(c_ranges...,n[1]-nh[1]+1:n[1]))
end

# Perform the copy without additional allocation
@generated function copy_ranges!{N}(c, d, c_ranges::NTuple{N}, d_ranges::NTuple{N})
    quote
        @nloops $N i x->1:length(c_ranges[x]) begin
         (@nref $N d x->d_ranges[x][i_x]) = (@nref $N c x->c_ranges[x][i_x])
        end
    end
end


#apply!{G,N,T}(op::Extension, dest, src::TensorProductBasis{FourierBasisOdd{T},G,N,T}, coef_dest::Array, coef_src::Array) = 
#    reshape_N_to_L!(coef_src, coef_dest, size(coef_src), size(coef_dest))
#
#apply!{G,N,T}(op::Restriction, dest::TensorProductBasis{FourierBasisOdd{T},G,N,T}, src, coef_dest::Array, coef_src::Array) =
#    reshape_L_to_N!(coef_dest, coef_src, size(coef_dest), size(coef_src))
#
#apply!{T,N,G,H,ELT}(op::Extension, dest, src::TensorProductBasis{DiscreteGridSpace1d{G,ELT,T},H,N,T}, coef_dest::Array, coef_src::Array)=reshape_N_to_L!(coef_src, coef_dest, size(coef_src), size(coef_dest))
#
#apply!{T,N,G,H,ELT}(op::Restriction, dest::TensorProductBasis{DiscreteGridSpace1d{G,ELT,T},H,N,T}, src, coef_dest::Array, coef_src::Array)=reshape_L_to_N!(coef_dest, coef_src, size(coef_dest), size(coef_src))


# TODO: fix proper dispatch here!
apply!(op::Extension, dest::TensorProductSet, src::TensorProductSet, coef_dest, coef_src) = 
    apply_tensor!(op, dest, src, set(dest,1), set(src,1), coef_dest, coef_src)

apply!(op::Restriction, dest::TensorProductSet, src::TensorProductSet, coef_dest, coef_src) =
    apply_tensor!(op, dest, src, set(dest,1), set(src,1), coef_dest, coef_src)

apply_tensor!(op::Extension, dest::TensorProductSet, src::TensorProductSet, dest1::FourierBasis, src1::FourierBasisOdd, coef_dest, coef_src) = 
    reshape_N_to_L!(coef_src, coef_dest, size(coef_src), size(coef_dest))

apply_tensor!(op::Restriction, dest::TensorProductSet, src::TensorProductSet, dest1::FourierBasisOdd, src1::FourierBasis, coef_dest, coef_src) =
    reshape_L_to_N!(coef_dest, coef_src, size(coef_dest), size(coef_src))

apply_tensor!(op::Extension, dest::TensorProductSet, src::TensorProductSet, dest1::DiscreteGridSpace, src1::DiscreteGridSpace, coef_dest, coef_src) = 
    reshape_N_to_L!(coef_src, coef_dest, size(coef_src), size(coef_dest))

apply_tensor!(op::Restriction, dest::TensorProductSet, src::TensorProductSet, dest1::DiscreteGridSpace, src1::DiscreteGridSpace, coef_dest, coef_src) =
    reshape_L_to_N!(coef_dest, coef_src, size(coef_dest), size(coef_src))

function fourier_extension_problem{T}(n::Int, t::T, sampling, domain::AbstractDomain1d{T})
    m = 2*round(Int, (n-1)/2 * sampling)+1
    l = round(Int, t*(m-1))
    fourier_extension_problem(n, m, l, domain)
end

function fourier_extension_problem{T}(n::Int, m::Int, l::Int, domain::Interval{T})
    @assert isodd(n)

    t = (l*one(T)) / ((m-1)*one(T))

    a = left(domain)
    b = right(domain)

    fbasis1 = FourierBasis(n, a, b + (b-a)*(t-1))
    fbasis2 = FourierBasis(l, a, b + (b-a)*(t-1))

    grid1 = grid(fbasis1)
    grid2 = grid(fbasis2)

    rgrid = EquispacedSubGrid(grid2, 1, m)

    tbasis1 = DiscreteGridSpace(grid1)
    tbasis2 = DiscreteGridSpace(grid2)

    tbasis_restricted = DiscreteGridSpace(rgrid)

    f_extension = Extension(fbasis1, fbasis2)
    f_restriction = Restriction(fbasis2, fbasis1)

    t_extension = Extension(tbasis_restricted, tbasis2)
    t_restriction = Restriction(tbasis2, tbasis_restricted)

    transform1 = transform_operator(tbasis1, fbasis1)
    itransform1 = transform_operator(fbasis1, tbasis1)

    transform2 = transform_operator(tbasis2, fbasis2)
    itransform2 = transform_operator(fbasis2, tbasis2)

    FE_DiscreteProblem(domain, fbasis1, fbasis2, tbasis1, tbasis2,
        tbasis_restricted, f_extension, f_restriction, t_extension,
        t_restriction, transform1, itransform1, transform2, itransform2)
end

function fourier_extension_problem{N,T}(n::NTuple{N,Int}, m::NTuple{N,Int}, l::NTuple{N,Int}, domain::Cube{N,T})   

    bbox=box(domain)
    fbasis1=Array{FourierBasisOdd}(N)
    fbasis2=Array{FourierBasisEven}(N)
    for i=1:N
        t = (l[1]*one(T))/((m[1]-1)*one(T))
        fbasis1[i] = FourierBasis(n[i], left(bbox)[i], right(bbox)[i] + (right(bbox)[i]-left(bbox)[i])*(t-1))
        fbasis2[i] = FourierBasis(l[i], left(bbox)[i], right(bbox)[i] + (right(bbox)[i]-left(bbox)[i])*(t-1))
    end
    tens_fbasis1 = TensorProductSet(ntuple(i->fbasis1[i], N)...)
    tens_fbasis2 = TensorProductSet(ntuple(i->fbasis2[i], N)...)

    tens_grid1 = TensorProductGrid(ntuple(i->grid(fbasis1[i]),N)...)
    tens_grid2 = TensorProductGrid(ntuple(i->grid(fbasis1[i]),N)...)

    tens_rgrid = TensorProductGrid(ntuple(i->EquispacedSubGrid(grid(fbasis2[i]), 1, m[i]), N)...)
    tens_tbasis1 = TensorProductSet(ntuple(i->DiscreteGridSpace(grid(fbasis1[i])), N)...)
    tens_tbasis2 = TensorProductSet(ntuple(i->DiscreteGridSpace(grid(fbasis2[i])), N)...)

    tens_tbasis_restricted = TensorProductSet(ntuple(i->DiscreteGridSpace(EquispacedSubGrid(grid(fbasis2[i]), 1, m[i])), N)...)

    f_extension = Extension(tens_fbasis1, tens_fbasis2)
    f_restriction = Restriction(tens_fbasis2, tens_fbasis1)

    t_extension = Extension(tens_tbasis_restricted, tens_tbasis2)
    t_restriction = Restriction(tens_tbasis2, tens_tbasis_restricted)

    transform1 = FastFourierTransformFFTW(tens_tbasis1, tens_fbasis1)
    itransform1 = InverseFastFourierTransformFFTW(tens_fbasis1, tens_tbasis1)

    transform2 = FastFourierTransformFFTW(tens_tbasis2, tens_fbasis2)
    itransform2 = InverseFastFourierTransformFFTW(tens_fbasis2, tens_tbasis2)

    FE_DiscreteProblem(domain, tens_fbasis1, tens_fbasis2, tens_tbasis1, tens_tbasis2,
        tens_tbasis_restricted, f_extension, f_restriction, t_extension,
        t_restriction, transform1, itransform1, transform2, itransform2)
end


function apply!{T,G <: MaskedGrid}(op::Extension, dest, src::DiscreteGridSpace{G}, coef_dest::Array{T}, coef_src::Array{T})
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)

    grid1 = grid(src)
    # Again too much work, but better than not filling at all
    fill!(coef_dest, zero(T))

    l = 0
    for i in eachindex(grid1.grid)
        if in(i, grid1)
            l = l+1
            coef_dest[i] = coef_src[l]
        end
    end

#    l = 0
#    for i in eachindex(grid1)
#        l = l+1
#        coef_dest[i] = coef_src[l]
#    end
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

#    l = 0
#    for i in eachindex(grid1)
#        l = l+1
#        coef_dest[l] = coef_src[i]
#    end
end


function fourier_extension_problem{N,T}(n::NTuple{N,Int}, m::NTuple{N,Int}, l::NTuple{N,Int},
                                        domain::AbstractDomain{N,T})

    le=left(box(domain))
    ri=right(box(domain))
    CubeBox = Cube(ntuple(i->le[i],N),ntuple(i->ri[i],N))
    problem = fourier_extension_problem(n, m, l, CubeBox)

    tbasis2 = problem.tbasis2
    tbasis_restricted = DiscreteGridSpace(MaskedGrid(grid(tbasis2), domain))

    t_extension = Extension(tbasis_restricted, tbasis2)
    t_restriction = Restriction(tbasis2, tbasis_restricted)

    FE_DiscreteProblem(domain, problem.fbasis1, problem.fbasis2, problem.tbasis1, problem.tbasis2, 
        tbasis_restricted, problem.f_extension, problem.f_restriction,
        t_extension, t_restriction,
        problem.transform1, problem.itransform1, problem.transform2, problem.itransform2)
end



######################
# Default parameters
######################

default_fourier_n(domain::AbstractDomain1d) = 50

default_fourier_n(domain::AbstractDomain2d) = (10, 10)

default_fourier_n(domain::AbstractDomain3d) = (3, 3, 3)

default_fourier_T{T}(domain::AbstractDomain{1,T}) = 2*one(T)
default_fourier_T{N,T}(domain::AbstractDomain{N,T}) = ntuple(i->2*one(T),N)

default_fourier_sampling{T}(domain::AbstractDomain{1,T}) = 2*one(T)
default_fourier_sampling{N,T}(domain::AbstractDomain{N,T}) = ntuple(i->2*one(T),N)


default_fourier_problem{T}(domain::AbstractDomain1d{T}, n::Int, t::Number, s::Number) =
    fourier_extension_problem(2*n+1, convert(T, t), convert(T,s), domain)

function default_fourier_problem{T}(domain::AbstractDomain2d{T}, n::Int, t::Number, s::Number)
    N = 2*n+1
    M = 2*round(Int, n*s)+1
    t = round(Int,t*(M-1)/2)*2/(M-1)
    L = round(Int, t*(M-1))
    fourier_extension_problem((N,N), (M,M), (L,L), domain)
end

function default_fourier_problem{T}(domain::AbstractDomain3d{T}, n::Int, t::Number, s::Number)
    N = 2*n+1
    M = 2*round(Int, n*s)+1
    t = round(Int,t*(M-1)/2)*2/(M-1)
    L = round(Int, t*(M-1))
    fourier_extension_problem((N,N,N), (M,M,M), (L,L,L), domain)
end

function default_fourier_problem{T}(domain::AbstractDomain{1,T}, n::Tuple, t::Tuple, s::Tuple)
    N = 2*[n...]+1
    M = 2*round(Int, [n...].*[s...])+1
    t = round(Int,[t...].*(M-1)/2).*(2./(M-1))
    L = round(Int, t.*(M-1))
    fourier_extension_problem(N[1],M[1],L[1], domain)
end

function default_fourier_problem{T}(domain::AbstractDomain{T}, n::Tuple, t::Tuple, s::Tuple)
    N = 2*[n...]+1
    M = 2*round(Int, [n...].*[s...])+1
    t = round(Int,[t...].*(M-1)/2).*(2./(M-1))
    L = round(Int, t.*(M-1))
    fourier_extension_problem(tuple(N...),tuple(M...),tuple(L...), domain)
end
#dirty

default_fourier_domain_1d() = Interval()

default_fourier_domain_2d() = Circle()

default_fourier_domain_3d() = Sphere()

default_fourier_solver(domain) = FE_DirectSolver

#default_fourier_solver(domain::Interval{Float64}) = FE_ProjectionSolver
default_fourier_solver(domain::Interval) = FE_ProjectionSolver

default_fourier_solver(domain::Cube) =  FE_TensorProductSolver
    


