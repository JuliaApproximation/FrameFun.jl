
immutable PlungeOperator{T,SRC,DEST} <: AbstractOperator{SRC,DEST}
    A       ::  AbstractOperator
    AstarA  ::  AbstractOperator
    lambda  ::  T

    PlungeOperator(op::FE_DiscreteOperator) = new(op, op*ctranspose(op), length(frequency_basis_ext(problem(op)))*one(numtype(op)))
end

PlungeOperator(op::FE_DiscreteOperator) = PlungeOperator{numtype(op),typeof(dest(op)),typeof(dest(op))}(op)

src(op::PlungeOperator) = src(op.AstarA)

dest(op::PlungeOperator) = dest(op.AstarA)

eltype(op::PlungeOperator) = promote_type(promote_type(eltype(src(op)),eltype(dest(op))),eltype(op.A))


function apply!(op::PlungeOperator, dest, src, coef_dest, coef_src)
    apply!(op.AstarA, coef_dest, coef_src)
    lambda = length(frequency_basis_ext(problem(op.A)))
    for i in eachindex(coef_dest)
        coef_dest[i] = coef_dest[i] - op.lambda*coef_src[i]
    end
end



immutable FE_ProjectionSolver <: FE_Solver
    op      ::  FE_DiscreteOperator
end


