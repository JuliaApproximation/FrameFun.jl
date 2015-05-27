# FrameFuns.jl

module FrameFuns

using BasisFunctions

using Base.Cartesian

import Base: length, eltype, size, push!

import Base: eachindex, start, next, done, getindex

import Base: show, showcompact

import BasisFunctions: src, dest, matrix, matrix!, apply!, numtype

import BasisFunctions: dim, index_dim, grid

import BasisFunctions: operator, coefficients, set


# from domains.jl
export Interval, Circle, Square, Cube, Rectangle, Cylinder, atomium

export numtype

# from funs.jl
export ExpFun


include("domains.jl")

include("subgrid.jl")

include("fe_problem.jl")

include("my_lsqr.jl")

include("fe_solvers.jl")

include("fastsolver.jl")

include("funs.jl")

include("fe_fourier.jl")

include("plots.jl")

end # module


