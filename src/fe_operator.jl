# fe_operator.jl

abstract FE_Operator{SRC,DEST} <: AbstractOperator{SRC,DEST}

# Delegation methods
for op in (:numtype,:eltype,:frequency_basis,:frequency_basis_ext)
    @eval $op(o::FE_Operator) = $op(problem(o))
end


problem{OP <: FE_Operator}(op::OperatorTranspose{OP}) = problem(operator(op))

# Delegation methods
for op in (:numtype,:eltype,:frequency_basis,:frequency_basis_ext)
    @eval $op{OP <: FE_Operator}(o::OperatorTranspose{OP}) = $op(problem(o))
end



immutable FE_DiscreteOperator{SRC,DEST} <: FE_Operator{SRC,DEST}
    problem ::  FE_DiscreteProblem

    scratch1
    scratch2        # For transforms that can be done in-place we don't need scratch2

    function FE_DiscreteOperator(problem)
        scratch1 = Array(eltype(problem), size(frequency_basis_ext(problem)))
        scratch2 = Array(eltype(problem), size(frequency_basis_ext(problem)))
        new(problem, scratch1, scratch2)
    end
end

FE_DiscreteOperator(problem::FE_DiscreteProblem) = FE_DiscreteOperator{typeof(frequency_basis(problem)),typeof(restricted_time_basis(problem))}(problem)

src(op::FE_DiscreteOperator) = frequency_basis(op.problem)

dest(op::FE_DiscreteOperator) = restricted_time_basis(op.problem)

problem(op::FE_DiscreteOperator) = op.problem




function apply!(op::FE_DiscreteOperator, dest, src, coef_dest, coef_src)
    p = problem(op)
#    apply!(p.f_extension, op.scratch1, coef_src)
#    apply!(p.itransform2, op.scratch2, op.scratch1)
#    apply!(p.t_restriction, coef_dest, op.scratch2)
    apply!(p.f_extension, op.scratch1, coef_src)
    apply!(p.itransform2, op.scratch1)
    apply!(p.t_restriction, coef_dest, op.scratch1)
end


function apply!(::OperatorTranspose, op::FE_DiscreteOperator, coef_dest, coef_src)
    p = problem(op)
#    apply!(p.t_extension, op.scratch2, coef_src)
#    apply!(p.transform2, op.scratch1, op.scratch2)
#    apply!(p.f_restriction, coef_dest, op.scratch1)
    apply!(p.t_extension, op.scratch1, coef_src)
    apply!(p.transform2, op.scratch1)
    apply!(p.f_restriction, coef_dest, op.scratch1)
end




