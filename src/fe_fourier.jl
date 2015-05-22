# fe_fourier.jl


function apply!{T}(op::ZeroPadding, dest::FourierBasis, src::FourierBasis, coef_dest::Array{T}, coef_src::Array{T})
    @assert size(coef_src)==size(src)
    @assert size(coef_dest)==size(dest)

    n = length(src)
    l = length(dest)

    nh = (n-1) >> 1

    for i = 1:nh+1
        coef_dest[i] = coef_src[i]
    end
    for i = nh+2:l-nh
        coef_dest[i] = zero(T)
    end
    for i = 1:nh
        coef_dest[end-nh+i] = coef_src[end-nh+i]
    end
end


function apply!(op::Restriction, dest::FourierBasis, src::FourierBasis, coef_dest::Array, coef_src::Array)
    @assert size(coef_src)==size(src)
    @assert size(coef_dest)==size(dest)

    n = length(dest)
    l = length(src)

    nh = (n-1) >> 1

    for i = 1:nh+1
        coef_dest[i] = coef_src[i]
    end
    for i = 1:nh
        coef_dest[end-nh+i] = coef_src[end-nh+i]
    end
end




# Reshape functions: we want to efficiently copy the data from a vector of length N to a larger vector of length L.
# Hard to do with Cartesian, but it can be done for any dimension recursively:
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


apply!{G,N,T}(op::ZeroPadding, dest, src::TensorProductBasis{FourierBasisOdd{T},G,N,T}, coef_dest::Array, coef_src::Array) = 
    reshape_N_to_L!(coef_src, coef_dest, size(coef_src), size(coef_dest))

apply!{G,N,T}(op::Restriction, dest::TensorProductBasis{FourierBasisOdd{T},G,N,T}, src, coef_dest::Array, coef_src::Array) =
    reshape_L_to_N!(coef_dest, coef_src, size(coef_dest), size(coef_src))

fourier_extension_problem{T <: FloatingPoint, S <: Number}(n::Int, t::T, sampling::S) = fourier_extension_problem(n, promote(t,sampling)...)

function fourier_extension_problem{T <: FloatingPoint}(n::Int, t::T, sampling::T)
    m = Int(ceil(n * sampling))
    l = Int(round(t*(m-1)))
    fourier_extension_problem(n, m, l)
end

function fourier_extension_problem{T}(n::Int, m::Int, l::Int, a::T = 0.0, b::T = 1.0)
    assert(isodd(n))

    t = (l*one(T)) / ((m-1)*one(T))

    domain = Interval(a, b)

    fbasis1 = FourierBasis(n, a, b + (b-a)*(t-1))
    fbasis2 = FourierBasis(l, a, b + (b-a)*(t-1))

    grid1 = grid(fbasis1)
    grid2 = grid(fbasis2)

    rgrid = EquispacedSubGrid(grid2, 1, m)

    tbasis1 = TimeDomain(grid1)
    tbasis2 = TimeDomain(grid2)

    tbasis_restricted = TimeDomain(rgrid)

    f_extension = ZeroPadding(fbasis1, fbasis2)
    f_restriction = Restriction(fbasis2, fbasis1)

    t_extension = ZeroPadding(tbasis_restricted, tbasis2)
    t_restriction = Restriction(tbasis2, tbasis_restricted)

    transform1 = FastFourierTransform(tbasis1, fbasis1)
    itransform1 = InverseFastFourierTransform(fbasis1, tbasis1)

    transform2 = FastFourierTransform(tbasis2, fbasis2)
    itransform2 = InverseFastFourierTransform(fbasis2, tbasis2)

    FE_DiscreteProblem(domain, fbasis1, fbasis2, tbasis1, tbasis2, tbasis_restricted, f_extension, f_restriction, t_extension, t_restriction, transform1, itransform1, transform2, itransform2)
end

function fourier_extension_problem{N}(n::NTuple{N,Int}, m::NTuple{N,Int}, l::NTuple{N,Int}, T = Float64)
#     T = BigFloat
#    T = Float64

    t = (l[1]*one(T)) / ((m[1]-1)*one(T))

    domain = Cube(N)

    fbasis1 = FourierBasis(n[1], -one(T), one(T) + 2*(t-1))
    fbasis2 = FourierBasis(l[1], -one(T), one(T) + 2*(t-1))
    tens_fbasis1 = tensorproduct(fbasis1, N)
    tens_fbasis2 = tensorproduct(fbasis2, N)

    grid1 = grid(fbasis1)
    grid2 = grid(fbasis2)
    tens_grid1 = tensorproduct(grid1, N)
    tens_grid2 = tensorproduct(grid2, N)

    rgrid = EquispacedSubGrid(grid2, 1, m[1])
    tens_rgrid = tensorproduct(rgrid, N)

    tbasis1 = TimeDomain(grid1)
    tbasis2 = TimeDomain(grid2)
    tens_tbasis1 = tensorproduct(tbasis1, N)
    tens_tbasis2 = tensorproduct(tbasis2, N)

    tbasis_restricted = TimeDomain(rgrid)
    tens_tbasis_restricted = tensorproduct(tbasis_restricted, N)

    f_extension = ZeroPadding(tens_fbasis1, tens_fbasis2)
    f_restriction = Restriction(tens_fbasis2, tens_fbasis1)

    t_extension = ZeroPadding(tens_tbasis_restricted, tens_tbasis2)
    t_restriction = Restriction(tens_tbasis2, tens_tbasis_restricted)

    transform1 = FastFourierTransform(tens_tbasis1, tens_fbasis1)
    itransform1 = InverseFastFourierTransform(tens_fbasis1, tens_tbasis1)

    transform2 = FastFourierTransform(tens_tbasis2, tens_fbasis2)
    itransform2 = InverseFastFourierTransform(tens_fbasis2, tens_tbasis2)

    FE_DiscreteProblem(domain, tens_fbasis1, tens_fbasis2, tens_tbasis1, tens_tbasis2,
        tens_tbasis_restricted, f_extension, f_restriction, t_extension,
        t_restriction, transform1, itransform1, transform2, itransform2)
end


function apply!{T,G <: MaskedGrid}(op::ZeroPadding, dest, src::TimeDomain{G}, coef_dest::Array{T}, coef_src::Array{T})
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)

    grid1 = grid(src)

    l = 0
    for i in eachindex_mask(grid1)
        if in(i, grid1)
            l = l+1
            coef_dest[i] = coef_src[l]
        else
            coef_dest[i] = 0
        end
    end
end


function apply!{T,G <: MaskedGrid}(op::Restriction, dest::TimeDomain{G}, src, coef_dest::Array{T}, coef_src::Array{T})
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)

    grid1 = grid(dest)

    l = 0
    for i in eachindex_mask(grid1)
        if in(i, grid1)
            l = l+1
            coef_dest[l] = coef_src[i]
        end
    end
end

function fourier_extension_problem{N,T}(n::NTuple{N,Int}, m::NTuple{N,Int}, l::NTuple{N,Int},
                                        domain::AbstractDomain{N,T})

    problem = fourier_extension_problem(n, m, l, T)
    
    tbasis2 = problem.tbasis2
    tbasis_restricted = TimeDomain(MaskedGrid(grid(tbasis2), domain))

    t_extension = ZeroPadding(tbasis_restricted, tbasis2)
    t_restriction = Restriction(tbasis2, tbasis_restricted)

    FE_DiscreteProblem(domain, problem.fbasis1, problem.fbasis2, problem.tbasis1, problem.tbasis2, 
        tbasis_restricted, problem.f_extension, problem.f_restriction,
        t_extension, t_restriction,
        problem.transform1, problem.itransform1, problem.transform2, problem.itransform2)
end


default_fourier_problem(domain::AbstractDomain1d) = fourier_extension_problem(21, 2.0, 2.0)

default_fourier_problem(domain::AbstractDomain2d) = fourier_extension_problem((11,11), (22,22), (44,44))

default_fourier_problem(domain::AbstractDomain3d) = fourier_extension_problem((5,5,5), (10,10,10), (20,20,20))

default_fourier_domain_1d() = Interval()

default_fourier_domain_2d() = Circle()

default_fourier_domain_3d() = Sphere()

default_fourier_solver(problem) = default_fourier_solver(problem, domain(problem))

default_fourier_solver(problem, domain::AbstractDomain) = FE_DirectSolver(problem)

default_fourier_solver(problem, domain::Interval) = FE_ProjectionSolver(problem)


