# fe_problem.jl

"""
An FE_Problem groups all the information of a Fourier extension problem.
This data can be used in a solver to produce an approximation.
"""
abstract FE_Problem{N,T}

numtype{N,T}(::Type{FE_Problem{N,T}}) = T
numtype{P <: FE_Problem}(::Type{P}) = numtype(super(P))
numtype(p::FE_Problem) = numtype(typeof(p))

dim{N,T}(::Type{FE_Problem{N,T}}) = N
dim{P <: FE_Problem}(::Type{P}) = dim(super(P))
dim(p::FE_Problem) = dim(typeof(p))

eltype(p::FE_Problem) = eltype(operator(p))



# This type groups the data corresponding to a FE problem.
immutable FE_DiscreteProblem{N,T} <: FE_Problem{N,T}
    domain          ::  AbstractDomain{N,T}

    op              ::  AbstractOperator
    opt             ::  AbstractOperator

    # The original and extended frequency basis
    fbasis1
    fbasis2
    # The original and extended time basis
    tbasis1
    tbasis2
    # Time domain basis restricted to the domain of the problem
    tbasis_restricted

    # Extension and restriction operators
    f_extension         # from fbasis1 to fbasis2
    f_restriction       # from fbasis2 to fbasis1
    t_extension         # from tbasis_restricted to tbasis2
    t_restriction       # from tbasis2 to tbasis_restricted

    # Transforms from time to frequency domain and back
    transform1          # from tbasis1 to fbasis1
    itransform1         # from fbasis1 to tbasis1
    transform2          # from tbasis2 to fbasis2
    itransform2         # from fbasis2 to tbasis2

    normalization       # transform normalization operator
end

#FE_DiscreteProblem{N,T}(domain::AbstractDomain{N,T}, otherargs...) =
#    FE_DiscreteProblem{N,T}(domain, otherargs...)


function FE_DiscreteProblem(domain::AbstractDomain, fbasis1, fbasis2, tbasis1, tbasis2, tbasis_restricted)
    ELT = promote_type(eltype(tbasis_restricted), eltype(fbasis1))
    f_extension = extension_operator(fbasis1, fbasis2, ELT)
    f_restriction = restriction_operator(fbasis2, fbasis1, ELT)

    t_extension = extension_operator(tbasis_restricted, tbasis2)
    t_restriction = restriction_operator(tbasis2, tbasis_restricted)

    transform1 = transform_operator(tbasis1, fbasis1)
    itransform1 = transform_operator(fbasis1, tbasis1)
    transform2 = transform_operator(tbasis2, fbasis2)
    itransform2 = transform_operator(fbasis2, tbasis2)

    normalization = transform_normalization_operator(fbasis1, fbasis2, ELT)

    op  = t_restriction * itransform2 * f_extension
    opt = f_restriction * transform2 * t_extension

    FE_DiscreteProblem(domain, op, opt, fbasis1, fbasis2, tbasis1, tbasis2, tbasis_restricted,
        f_extension, f_restriction, t_extension, t_restriction,
        transform1, itransform1, transform2, itransform2, normalization)
end



domain(p::FE_DiscreteProblem) = p.domain

operator(p::FE_DiscreteProblem) = p.op

operator_transpose(p::FE_DiscreteProblem) = p.opt

normalization(p::FE_DiscreteProblem) = p.normalization

frequency_basis(p::FE_DiscreteProblem) = p.fbasis1
frequency_basis_ext(p::FE_DiscreteProblem) = p.fbasis2

time_basis(p::FE_DiscreteProblem) = p.tbasis1
time_basis_ext(p::FE_DiscreteProblem) = p.tbasis2
time_basis_restricted(p::FE_DiscreteProblem) = p.tbasis_restricted


f_extension(p::FE_DiscreteProblem) = p.f_extension
f_restriction(p::FE_DiscreteProblem) = p.f_restriction

t_extension(p::FE_DiscreteProblem) = p.t_extension
t_restriction(p::FE_DiscreteProblem) = p.t_restriction


transform1(p::FE_DiscreteProblem) = p.transform1
itransform1(p::FE_DiscreteProblem) = p.itransform1
transform2(p::FE_DiscreteProblem) = p.transform2
itransform2(p::FE_DiscreteProblem) = p.itransform2


size(p::FE_DiscreteProblem) = size(operator(p))

size(p::FE_DiscreteProblem, j) = size(operator(p), j)

size_ext(p::FE_DiscreteProblem) = size(frequency_basis_ext(p))

length_ext(p::FE_DiscreteProblem) = length(frequency_basis_ext(p))

size_ext(p::FE_DiscreteProblem, j) = size(frequency_basis_ext(p), j)

param_N(p::FE_DiscreteProblem) = length(frequency_basis(p))

param_L(p::FE_DiscreteProblem) = length(time_basis_ext(p))

param_M(p::FE_DiscreteProblem) = length(time_basis_restricted(p))



function rhs(p::FE_Problem, f::Function, elt = eltype(p))
    grid1 = grid(time_basis_restricted(p))
    M = length(grid1)
    b = Array(elt, M)
    rhs!(p, b, f)
    b
end

function rhs!(p::FE_Problem, b::AbstractArray, f::Function)
    grid1 = grid(time_basis_restricted(p))
    M = length(grid1)

    @assert length(b) == M

    rhs!(grid1, b, f)
end

function rhs!(grid::AbstractGrid, b::AbstractArray, f::Function)
    for i in eachindex(grid)
        b[i] = f(grid[i]...)
    end
end


immutable FE_TensorProductProblem{TP,PN,N,T} <: FE_Problem{N,T}
    problems       ::  TP
end

function FE_TensorProductProblem(problems...)
        TP = typeof(problems)
        PN = length(problems)
        T = numtype(problems[1])
        N = sum(map(dim, problems))
        FE_TensorProductProblem{TP,PN,N,T}(problems)
end


domain(p::FE_TensorProductProblem) = TensorProductDomain(map(domain,p.problems)...)

tp_length{TP,PN,N,T}(p::FE_TensorProductProblem{TP,PN,N,T}) = PN

for op in (:frequency_basis, :frequency_basis_ext, :time_basis_restricted, :time_basis_ext)
    @eval $op{TP,PN}(p::FE_TensorProductProblem{TP,PN}) = 
        TensorProductSet(map($op,p.problems)...)
end

for op in (:operator, :operator_transpose, :t_restriction, :t_extension,
           :itransform2, :transform2, :f_restriction, :f_extension,
           :normalization_operator)
    @eval $op{TP,PN}(p::FE_TensorProductProblem{TP,PN}) =
        TensorProductOperator(map($op,p.problems)...)
end

for op in (:dim, :length_ext, :param_N, :param_L, :param_M)
    @eval $op{TP,PN}(p::FE_TensorProductProblem{TP,PN}) = 
        prod(map($op,p.problems)...)
end

for op in (:size, :size_ext)
    @eval $op{TP,PN}(p::FE_TensorProductProblem{TP,PN}) = 
        tuple(map($op,p.problems)...)
end


# This code needs a revision.
# What is the true (generic) meaning of transform_normalization_operator when there is a source and a destination?
# We also have to do something about the types added as arguments here.
transform_normalization_operator(src::TensorProductSet, dest::TensorProductSet, ELT) =
    TensorProductOperator(ELT, [transform_normalization_operator(set(src,i), set(dest,i), ELT) for i in 1:tp_length(src)]...)

function transform_normalization_operator(src::FunctionSet, dest::FunctionSet, ELT)
    factor = sqrt(ELT(length(src))/ELT(length(dest)))
    transform_normalization_operator(src, ELT) * ScalingOperator(src, factor)
end

# Perhaps this is not always correct. Check.
function transform_normalization_operator(p::FE_DiscreteProblem, ELT)
    transform_normalization_operator(frequency_basis(p), frequency_basis_ext(p), ELT)
end

function transform_normalization_operator(p::FE_TensorProductProblem, ELT)
    TensorProductOperator(ELT, [transform_normalization_operator(p.problems[i], ELT) for i in 1:tp_length(p)]...)
end



