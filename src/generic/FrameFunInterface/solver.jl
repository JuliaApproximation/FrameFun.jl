include("solvers/azsolver.jl")
include("solvers/generic_azsolver.jl")
include("solvers/lowranksolver.jl")
include("solvers/randomized.jl")
include("solvers/smoothsolver.jl")

# solver needs also a problemstyle and a solverstyle
solver(ap::ApproximationProblem;
            problemstyle = ProblemStyle(ap),
            samplingstyle = SamplingStyle(ap),
            solverstyle = SolverStyle(samplingstyle, ap), options...) =
        solver(problemstyle, solverstyle, ap; samplingstyle=samplingstyle, options...)




# solver needs a solverstyle ap and operator
function solver(pstyle::DictionaryOperatorStyle, solverstyle::SolverStyle, ap::ApproximationProblem;
            samplingstyle = SamplingStyle(ap), options...)
    S = samplingoperator(samplingstyle, ap; options...)
    A = discretization(samplingstyle, ap, S; options...)
    solver(solverstyle, ap, A; S = S, options...)
end

# samplingoperator is no TensorProductOperator
function solver(pstyle::DictionaryOperatorStyle, solverstyle::ProductSolverStyle, ap::ApproximationProblem;
            samplingstyle = SamplingStyle(ap), options...)
    S = samplingoperator(samplingstyle, ap; options...)
    A = discretization(samplingstyle, ap; options...)
    solver(solverstyle, ap, A; S = S, options...)
end

solver(::InverseStyle, ap::ApproximationProblem, A::AbstractOperator; options...) = inv(A)
solver(::DirectStyle, ap::ApproximationProblem, A::AbstractOperator; options...) = directsolver(A; options...)
solver(::IterativeStyle, ap::ApproximationProblem, A::AbstractOperator; options...) = iterativesolver(A; options...)






solver(style::AZStyle, ap::ApproximationProblem, A::AbstractOperator; problemstyle=ProblemStyle(ap), options...) =
    solver(problemstyle, style, ap, A; options...)
solver(pstyle::GenericOperatorStyle, solverstyle::AZStyle, ap::ApproximationProblem, A; options...) =
    GenericAZSolver(ap, A; options...)
solver(pstyle::DictionaryOperatorStyle, style::AZStyle, ap::ApproximationProblem, A::AbstractOperator; options...) =
    solver(style, ap, A, AZ_Zt(pstyle, ap; (options)...); options...)



function solver(::AZStyle, ap::ApproximationProblem, A::AbstractOperator, Zt::AbstractOperator;
            B=nothing, smallcoefficients=false, smallcoefficients_atol=NaN, smallcoefficients_rtol=NaN, verbose=false, options...)
    if smallcoefficients
        w = BasisFunctions.gaussweights(sampling_grid(ap; options...), measure(ap; options...))
        normF = abs(sqrt(sum(w .* B.^2)))
        if !isnan(smallcoefficients_rtol)
            verbose && println("Change smallcoefficients relative tolerance to absolute tolerance rtol*||f||")
            smallcoefficients_atol = smallcoefficients_rtol*normF
            smallcoefficients_rtol = NaN
        end

        AZSolver(A, Zt; smallcoefficients=smallcoefficients, smallcoefficients_rtol=smallcoefficients_rtol,
                smallcoefficients_atol=smallcoefficients_atol, verbose=verbose, options...)
    else
        AZSolver(A, Zt; B=B, verbose=verbose, options...)
    end
end

solver(style::AZSmoothStyle, ap::ApproximationProblem, A::AbstractOperator; options...) = solver(style, ap, A, AZ_Zt(ap; options...); options...)
solver(::AZSmoothStyle, ap::ApproximationProblem, A::AbstractOperator, Zt; options...) = AZSolver_with_smoothing(A, Zt; options...)

solver(::DualStyle, ap::ApproximationProblem, A::AbstractOperator; options...) = dualdiscretization(ap; options...)'

solver(solverstyle::ProductSolverStyle, ap::ApproximationProblem, A::AbstractOperator; samplingstyle=SamplingStyle(ap), options...) =
    solver(solverstyle, samplingstyle, ap, A; options...)
function solver(solverstyle::ProductSolverStyle, samplingstyle::ProductSamplingStyle, ap::ApproximationProblem, A::AbstractOperator; S, options...)
    solvere, sse, ape, Ae, Se = elements(solverstyle), elements(samplingstyle), elements(ap), elements(A), BasisFunctions.productelements(S)
    @assert length(solvere) == length(sse) == length(ape) == length(Ae) == length(Se)
    TensorProductOperator(
		map( (solversi, ssi, api, Ai, Si)->solver(solversi, api, Ai; S=Si, samplingstyle=ssi, options...), solvere, sse, ape, Ae, Se)...
	)
end
