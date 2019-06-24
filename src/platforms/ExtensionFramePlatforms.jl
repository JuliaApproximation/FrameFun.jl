module ExtensionFramePlatforms
using ..Platforms
import ..Platforms: platform, SolverStyle, SamplingStyle, measure, dictionary, dualdictionary
using ..ExtensionFrames
import ..ExtensionFrames: support
using DomainSets, BasisFunctions

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

support(platform::ExtensionFramePlatform) = platform.domain

SolverStyle(dict::ExtensionFrame, ::OversamplingStyle) = AZStyle()

SamplingStyle(p::ExtensionFramePlatform) = OversamplingStyle()

dictionary(p::ExtensionFramePlatform, n) =
    extensionframe(p.domain, dictionary(p.basisplatform, n))

measure(platform::ExtensionFramePlatform) = restrict(measure(platform.basisplatform), platform.domain)

dualdictionary(platform::ExtensionFramePlatform, param, measure::Measure; options...) =
   extensionframe(dualdictionary(platform.basisplatform, param, supermeasure(measure); options...), platform.domain)

end
