module FrameFun

using RecipesBase, FillArrays, Reexport
using DomainSets, GridArrays
@reexport using BasisFunctions
# Submodules
include("FrameFunDomains//FrameFunDomains.jl")
@reexport using .FrameFunDomains
include("generic/Platforms/Platforms.jl")
@reexport using .Platforms

include("ExtensionFrames.jl")
@reexport using .ExtensionFrames

include("generic/ApproximationProblems.jl")
using .ApproximationProblems

include("fun/DictFuns.jl")
@reexport using .DictFuns

include("generic/FrameFunInterface/FrameFunInterface.jl")
@reexport using .FrameFunInterface


# Other platforms
include("platforms/AugmentationPlatforms.jl")
include("platforms/ExtensionFramePlatforms.jl")
include("platforms/WeightedSumPlatforms.jl")
include("platforms/FrameFunPlatforms.jl")
@reexport using .AugmentationPlatforms, .ExtensionFramePlatforms, .WeightedSumPlatforms, .FrameFunPlatforms

include("FrameFunInterfaceExtension/FrameFunInterfaceExtension.jl")
using .FrameFunInterfaceExtension



export residual, abserror, maxerror, L2error
include("fun/error.jl")

include("applications/DiffEquations.jl")
@reexport using .DiffEquations
include("applications/WeightedApproximation.jl")
@reexport using .WeightedApproximation

include("generic/Adaptivity.jl")
@reexport using .Adaptivity

include("recipes.jl")

end # module
