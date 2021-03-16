module ExtensionFramePlatforms
using ..Platforms
import ..Platforms: platform, SolverStyle, SamplingStyle, measure, dictionary,
    dualdictionary, correctparamformat, unsafe_dictionary
using ..ExtensionFrames
import ..ExtensionFrames: support
using DomainSets, BasisFunctions
using BasisFunctions: Measure

platform(dict::ExtensionFrame) = ExtensionFramePlatform(platform(basis(dict)), support(dict))

export ExtensionFramePlatform
"""
    struct ExtensionFramePlatform <: FramePlatform

A platform that creates extension frames.
"""
struct ExtensionFramePlatform <: FramePlatform
    basisplatform   ::  Platform
    domain          ::  Domain
end

correctparamformat(platform::ExtensionFramePlatform, param) =
    correctparamformat(platform.basisplatform, param)

support(platform::ExtensionFramePlatform) = platform.domain

SolverStyle(dict::ExtensionFrame, ::OversamplingStyle) = AZStyle()

SamplingStyle(p::ExtensionFramePlatform) = OversamplingStyle()

unsafe_dictionary(p::ExtensionFramePlatform, n) =
    extensionframe(p.domain, unsafe_dictionary(p.basisplatform, n))

measure(platform::ExtensionFramePlatform) = restrict(measure(platform.basisplatform), platform.domain)

dualdictionary(platform::ExtensionFramePlatform, param, measure::Measure; options...) =
   extensionframe(dualdictionary(platform.basisplatform, param, supermeasure(measure); options...), platform.domain)


import ..ParameterPaths: default_param_path

default_param_path(platform::ExtensionFramePlatform) =
   default_param_path(platform.basisplatform)


end
