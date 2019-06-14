
"""
A `ProductPlatform` corresponds to the product of two or more platforms. It results in a product
dictionary, product sampling operator, product solver, ... etcetera.
"""
struct ProductPlatform{N} <: Platform
    platforms   ::  NTuple{N,Platform}
end

ProductPlatform(platforms::Platform...) = ProductPlatform(platforms)

ProductPlatform(platform::Platform, n::Int) = ProductPlatform(ntuple(x->platform, n)...)

iscomposite(p::ProductPlatform) = true
elements(p::ProductPlatform) = p.platforms

platform(p::ProductPlatform, i) = p.platforms[i]

productparameter(p::ProductPlatform{N}, n::Int) where {N} = ntuple(x->n, N)
productparameter(p::ProductPlatform{N}, n::NTuple{N,Any}) where {N} = n

SamplingStyle(platform::ProductPlatform) = ProductSamplingStyle(map(SamplingStyle, elements(platform)))
SolverStyle(p::ProductPlatform, samplingstyle) =
    ProductSolverStyle(map(SolverStyle, elements(p), elements(samplingstyle)))

dictionary(p::ProductPlatform, n::Int) = dictionary(p, productparameter(p, n))
azdual_dict(sstyle::DiscreteStyle, p::ProductPlatform, n::Int; options...) = azdual_dict(sstyle, p, productparameter(p, n); options...)

dictionary(p::ProductPlatform, n) = tensorproduct(map(dictionary, elements(p), n)...)
azdual_dict(sstyle::ProductSamplingStyle, p::ProductPlatform, n; options...) = tensorproduct(map((x,y,z)->azdual_dict(x,y,z; options...), elements(sstyle), elements(p), n)...)

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
mutable struct ProductPlatformApproximation{N} <: ApproximationProblem
    platform        ::  ProductPlatform{N}
    param
    productparam    ::  NTuple{N,Any}
    dict            ::  TensorProductDict
    samplingparam

    ProductPlatformApproximation{N}(platform, param) where {N} =
        new(platform, param, productparameter(platform, param), dictionary(platform, param))
end

approximationproblem(platform::ProductPlatform{N}, param) where {N} =
    ProductPlatformApproximation{N}(platform, param)

SamplingStyle(ap::ProductPlatformApproximation) = SamplingStyle(ap.platform)
SolverStyle(samplingstyle::SamplingStyle, ap::ProductPlatformApproximation) = SolverStyle(ap.platform, samplingstyle)

elements(ap::ProductPlatformApproximation) = map(approximationproblem, elements(ap.platform), ap.productparam)

dictionary(ap::ProductPlatformApproximation) = ap.dict
