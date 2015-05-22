# FrameFuns.jl

module FrameFuns

using BasisFunctions

using Base.Cartesian

import Base: length, eltype, size, in, eachindex, getindex, push!

import Base: show, showcompact

import BasisFunctions: src, dest, matrix, matrix!, apply!, numtype

import BasisFunctions: index_dim, grid

import BasisFunctions: operator, coefficients, set

import BasisFunctions: plot



# from domains.jl
export Interval, Circle, Square, Cube, Rectangle, Cylinder, atomium

# from funs.jl
export ExpFun


include("domains.jl")

include("funs.jl")

include("subgrid.jl")

include("fe_problem.jl")

include("my_lsqr.jl")

include("fe_solvers.jl")

include("fe_fourier.jl")

include("fastsolver.jl")

include("plots.jl")

end # module


