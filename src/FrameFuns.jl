# FrameFuns.jl

module FrameFuns

using BasisFunctions

using Base.Cartesian

import Base: length, eltype, size, in

import BasisFunctions: src, dest, operator_matrix, apply!, numtype

export ExpFun

include("domains.jl")

include("funs.jl")

include("subgrid.jl")

include("fe_problem.jl")

include("fe_operator.jl")

include("fe_solvers.jl")

include("fe_fourier.jl")


end # module


