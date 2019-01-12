
struct ProductPlatform{N} <: Platform
    platforms   ::  NTuple{N,Platform}
end

ProductPlatform(platforms::Platform...) = ProductPlatform(platforms)

iscomposite(p::ProductPlatform) = true
elements(p::ProductPlatform) = p.platforms

platform(p::ProductPlatform, i) = p.platforms[i]

productparameter(p::ProductPlatform{N}, n::Int) where {N} = ntuple(x->n, N)
productparameter(p::ProductPlatform{N}, n::NTuple{N,Any}) where {N} = n

SolverStyle(p::ProductPlatform, dstyle) = ProductSolverStyle()

dictionary(p::ProductPlatform, n::Int) = dictionary(p, productparameter(p, n))
dualdictionary(p::ProductPlatform, n::Int) = dualdictionary(p, productparameter(p, n))

dictionary(p::ProductPlatform, n) = tensorproduct(map(dictionary, elements(p), n)...)
dualdictionary(p::ProductPlatform, n) = tensorproduct(map(dualdictionary, elements(p), n)...)

discretization(p::ProductPlatform, S; options...) =
    tensorproduct(map(discretization, elements(p))...)
dualdiscretization(p::ProductPlatform, S; options...) =
    tensorproduct(map(discretization, elements(p))...)

samplingoperator(p::ProductPlatform, n::Int; options...) =
    samplingoperator(p, productparameter(p, n); options...)
samplingoperator(p::ProductPlatform, param; options...) =
    tensorproduct(map(samplingoperator, elements(p), param)...)

dualsamplingoperator(p::ProductPlatform, n::Int, m::Int; options...) =
    dualsamplingoperator(p, productparameter(p, n), productparameter(p, m); options...)
dualsamplingoperator(p::ProductPlatform, n, m; options...) =
    tensorproduct(map(dualsamplingoperator, elements(p), n, m)...)


##############################
# Approximation problem type
##############################

"""
A `ProductPlatformApproximation` corresponds to a product platform.
"""
struct ProductPlatformApproximation{N} <: ApproximationProblem
    platform        ::  ProductPlatform{N}
    param
    productparam    ::  NTuple{N,Any}
    dict            ::  TensorProductDict

    ProductPlatformApproximation{N}(platform, param) where {N} =
        new(platform, param, productparameter(platfor, param), dictionary(platform, param))
end

approximationproblem(platform::ProductPlatform{N}, param) where {N} =
    ProductPlatformApproximation{N}(platform, param)
