module FrameFun

using FillArrays, GridArrays, StaticArrays
using LinearAlgebra, LowRankApprox
using Printf, Reexport
using SpecialFunctions
using DomainSets
using QuadGK
using CompositeTypes, CompositeTypes.Display

import CompositeTypes: components, component, ncomponents
import DomainSets: factors, nfactors, factor

@reexport using BasisFunctions
import BasisFunctions: coefficienttype, dual, weight, expansion

include("extra/FrameFunDomains/FrameFunDomains.jl")
@reexport using .FrameFunDomains

include("platforms/generic/platforms.jl")
include("platforms/generic/parameterpaths.jl")

include("frames/extensionframes.jl")

include("generic/approximationproblems.jl")
include("generic/interface/interface.jl")

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
# include("applications/WeightedApproximation.jl")
# @reexport using .WeightedApproximation
include("applications/high_dimensional.jl")

include("generic/adaptivity.jl")
include("generic/funs.jl")

include("util/plot.jl")

export chebvariable, chebvariables,
    funvariable, funvariables

end # module
