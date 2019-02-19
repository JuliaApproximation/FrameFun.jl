struct GenericAZSolver <: AbstractSolverOperator
    A           ::  SynthesisOperator
    Zt          ::  AnalysisOperator
    plunge_op   ::  AbstractOperator     # The plunge operator P = (I-A*Zt)
    S           ::  SamplingOperator     # The sampling operator S
    psolver     ::  DictionaryOperator   # The low rank plunge solver for S*P*A
end

function GenericAZSolver(ap::ApproximationProblem, A;
            REG = default_regularization,
            rankestimate = 40,
            threshold = default_threshold(A),
            options...)
    pstyle = GenericOperatorStyle()
    Zt = AZ_Zt(pstyle, ap; options...)
    P = plungeoperator(pstyle, ap; options...)
    M = plungematrix(pstyle, ap; options...) # TODO some inner products are calculated twice.
    S = samplingoperator(ap; options...)

    psolver = REG(M; threshold = threshold, rankestimate = rankestimate, options...)

    GenericAZSolver(A, Zt, P, S, psolver)
end


default_threshold(A::SynthesisOperator) = regularization_threshold(coefficienttype(src(A)))

operator(op::GenericAZSolver) = op.A

function apply(op::GenericAZSolver, fun; options...)
    Pb = apply(op.S, apply(op.plunge_op, fun; options...); options...)
    x1 = op.psolver*Pb
    fun1 = x->fun(x...)-(op.A*x1)(x...)
    x2 = apply(op.Zt, fun1; options...)
    x1+x2
end
