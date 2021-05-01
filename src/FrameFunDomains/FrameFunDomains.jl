module FrameFunDomains

using DomainSets, GridArrays, StaticArrays

import Base: show, isapprox, in, <, <=, >, >=
import DomainSets: indomain
import GridArrays: boundingbox
include("extensions.jl")
include("fourierdomains.jl")
include("fractals.jl")
include("atomium.jl")
include("characteristic.jl")

include("polardomain.jl")
end
