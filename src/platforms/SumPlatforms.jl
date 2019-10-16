module SumPlatforms
using ..Platforms
using  BasisFunctions
using BasisFunctions: ZeroOperator, Domain
import ..Platforms: SolverStyle, dictionary, dualdictionary, measure, param_double,
    unsafe_dictionary, correctparamformat

import BasisFunctions: elements, size, length, support,
    grid_evaluation_operator

export SumPlatform
"""
    struct SumPlatform <: FramePlatform

An sum platform is the union of two platforms.
"""
struct SumPlatform <: FramePlatform
    platform1       ::  Platform
    platform2       ::  Platform

    param           ::  Function
    SumPlatform(plat1::Platform,plat2::Platform) =
        new(plat1,plat2,k->(k,k))

    SumPlatform(plat1::Platform,plat2::Platform, f) =
        new(plat1,plat2,f)
end

struct ZeroDict{S,T} <: Dictionary{S,T}
    n::Int
    domain::Domain
end

size(dict::ZeroDict) = (dict.n,)
length(dict::ZeroDict) = dict.n
support(dict::ZeroDict) = dict.domain
elements(platform::SumPlatform) = (platform.platform1,platform.platform2)
SolverStyle(platform::SumPlatform, ::SamplingStyle) = AZStyle()
correctparamformat(platform::SumPlatform, param) =
    correctparamformat(platform.platform1, param) && correctparamformat(platform.platform2, param)
unsafe_dictionary(platform::SumPlatform, param) =
    MultiDict([map(dictionary, elements(platform), platform.param(param))...])

function dualdictionary(platform::SumPlatform, param, measure::Measure; options...)
    param1, param2 = platform.param(param)
    dict1 = dualdictionary(platform.platform1, param1, measure)
    dict2 = dualdictionary(platform.platform2, param2, measure)

    # dict2 = ZeroDict{domaintype(dict1),coefficienttype(dict1)}(param2,support(dict1))
    multidict(dict1,dict2)
end

grid_evaluation_operator(dict::ZeroDict, gb::GridBasis, grid; T = op_eltype(dict, gb), options...) =
    ZeroOperator{T}(GridBasis(dict, grid), dict)

measure(platform::SumPlatform) = measure(platform.platform1)

end
