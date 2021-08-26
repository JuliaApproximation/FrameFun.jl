module FrameFun

using RecipesBase, FillArrays, Reexport
using DomainSets, GridArrays

@reexport using BasisFunctions

include("extra/FrameFunDomains/FrameFunDomains.jl")
@reexport using .FrameFunDomains

include("platforms/generic/platforms.jl")
include("platforms/generic/parameterpaths.jl")

include("frames/extensionframes.jl")

include("generic/approximationproblems.jl")
include("generic/interface/FrameFunInterface.jl")

# Other platforms
include("platforms/augmentation.jl")
include("platforms/extensionplatform.jl")
include("platforms/weightedsumplatform.jl")
include("platforms/sumplatform.jl")
include("platforms/basesandframes.jl")

include("generic/interface/interface_extension.jl")

export residual, abserror, maxerror, L2error
include("generic/error.jl")

include("applications/DiffEquations.jl")
@reexport using .DiffEquations
include("applications/WeightedApproximation.jl")
@reexport using .WeightedApproximation
include("applications/high_dimensional.jl")

include("generic/adaptivity.jl")

include("recipes.jl")

end # module
