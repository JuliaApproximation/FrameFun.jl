
using BasisFunctions, Test

using BasisFunctions: iscompatible

@testset "FrameFunDomains" begin
    include("test_framefundomains.jl")
end
@testset "Platforms" begin
    include("test_platforms.jl")
end
@testset "ExtensionFrames" begin
    include("test_extensionframes.jl")
end
@testset "FrameFunInterface" begin
    include("test_framefuninterface.jl")
end
@testset "ExtensionFramePlatforms" begin
    include("test_extensionframeplatforms.jl")
end
@testset "Other Frame Platforms" begin
    include("test_otherframeplatforms.jl")
end
@testset "FrameFunPlatforms" begin
    include("test_framefunplatforms.jl")
end
@testset "FrameFunInterfaceExtension" begin
    include("test_framefuninterfaceextension.jl")
end
@testset "Adaptivity" begin
    include("test_adaptivity.jl")
end
include("test_scenariolist.jl")
include("test_suite.jl")
include("test_suite_applications.jl")
include("test_specificplatforms.jl")
include("test_plots.jl")

# include("test_notebooks.jl")
# include("create_readme.jl")


println("All tests succeeded!!!")
