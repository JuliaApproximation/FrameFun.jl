# fe_operator.jl

abstract FE_Operator{SRC,DEST} <: AbstractOperator{SRC,DEST}

numtype(op::FE_Operator) = numtype(fe_problem(op))


frequency_basis(op::FE_Operator) = frequency_basis(fe_problem(op))


immutable FE_DiscreteOperator{SRC,DEST} <: FE_Operator{SRC,DEST}
    problem ::  FE_DiscreteProblem
end

FE_DiscreteOperator(problem::FE_DiscreteProblem) = FE_DiscreteOperator{typeof(frequency_basis(problem)),typeof(restricted_time_basis(problem))}(problem)

src(op::FE_DiscreteOperator) = frequency_basis(op.problem)

dest(op::FE_DiscreteOperator) = restricted_time_basis(op.problem)

fe_problem(op::FE_DiscreteOperator) = op.problem


function apply!(op::FE_DiscreteOperator, coef_dest, coef_src)
    p = op.problem
    apply!(p.f_extension, p.scratch1, coef_src)
    apply!(p.itransform2, p.scratch2, p.scratch1)
    apply!(p.t_restriction, coef_dest, p.scratch2)
end


fe_matrix(op::FE_DiscreteOperator) = operator_matrix(op)

fe_matrix!(op::FE_DiscreteOperator, matrix) = operator_matrix!(op, matrix)

function fe_rhs(op::FE_DiscreteOperator, f::Function, elt = eltype(op))
    grid = natural_grid(restricted_time_basis(fe_problem(op)))
    M = length(grid)
    b = Array(elt, M)
    fe_rhs!(op, b, f)
    b
end

function fe_rhs!(op::FE_DiscreteOperator, b::Vector, f::Function)
    grid = natural_grid(restricted_time_basis(fe_problem(op)))
    M = length(grid)

    @assert length(b) == M

    for i = 1:M
        b[i] = f(grid[i])
    end
end


function apply!{G <: MaskedGrid,T}(op::ZeroPadding, dest::TimeDomain, src::TimeDomain{G}, coef_dest::Array{T}, coef_src::Array{T})
    
end


