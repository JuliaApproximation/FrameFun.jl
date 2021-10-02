
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

dualdictionary(platform::ExtensionFramePlatform, plt_par, measure::Measure; options...) =
   extensionframe(dualdictionary(platform.basisplatform, plt_par, supermeasure(measure); options...), platform.domain)


default_param_path(platform::ExtensionFramePlatform) =
   default_param_path(platform.basisplatform)
