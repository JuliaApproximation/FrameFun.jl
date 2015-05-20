# FrameFuns.jl

module FrameFuns

using BasisFunctions

using Base.Cartesian

import Base: length, eltype, size, in, eachindex, getindex

import BasisFunctions: src, dest, matrix, matrix!, apply!, numtype

import BasisFunctions: index_dim, grid

import BasisFunctions: operator

export ExpFun

include("domains.jl")

include("funs.jl")

include("subgrid.jl")

include("fe_problem.jl")

include("fe_operator.jl")

include("fe_solvers.jl")

include("fe_fourier.jl")


end # module


