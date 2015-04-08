# FE.jl

module FE

using BasisFunctions

using Base.Cartesian

import Base: length, eltype, size, in

import BasisFunctions: src, dest, operator_matrix, apply!, numtype

include("domains.jl")

include("subgrid.jl")

include("fe_problem.jl")

include("fe_operator.jl")

include("fe_solvers.jl")

include("fe_fourier.jl")


end # module


