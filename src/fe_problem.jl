# fe_problem.jl

abstract FE_Problem

# This type groups the data corresponding to a FE problem.
immutable FE_DiscreteProblem <: FE_Problem
    fbasis1
    fbasis2

    tbasis1
    tbasis2

    restricted_tbasis

    f_extension
    f_restriction

    t_extension
    t_restriction

    transform1
    itransform1

    transform2
    itransform2

    scratch1
    scratch2
end


numtype(p::FE_DiscreteProblem) = numtype(p.fbasis1)

eltype(p::FE_DiscreteProblem) = eltype(frequency_basis(p))

dim(p::FE_DiscreteProblem) = dim(p.fbasis1)

frequency_basis(p::FE_DiscreteProblem) = p.fbasis1

frequency_basis_ext(p::FE_DiscreteProblem) = p.fbasis2

time_basis(p::FE_DiscreteProblem) = p.tbasis1

time_basis_ext(p::FE_DiscreteProblem) = p.tbasis2

restricted_time_basis(p::FE_DiscreteProblem) = p.restricted_tbasis




