include("../solvers/azsolver.jl")
include("../solvers/lowranksolver.jl")
include("../solvers/randomized.jl")
include("../solvers/smoothsolver.jl")

# Dispatch on solver style
solver(ap::ApproximationProblem; options...) = solver(SolverStyle(ap), ap; options...)

# add the linear system as a third argument
function solver(solverstyle::SolverStyle, ap::ApproximationProblem; options...)
    A = discretization(ap; options...)
    solver(solverstyle, ap, A; options...)
end

solver(::InverseStyle, ap::ApproximationProblem, A::AbstractOperator; options...) = inv(A)
solver(::DirectStyle, ap::ApproximationProblem, A::AbstractOperator; options...) =
	directsolver(A; options...)
solver(::IterativeStyle, ap::ApproximationProblem, A::AbstractOperator; options...) =
    iterativesolver(A; options...)
solver(::DualStyle, ap::ApproximationProblem, A::AbstractOperator; options...) =
	dualdiscretization(ap; options...)'

solver(solverstyle::ProductSolverStyle, ap::ApproximationProblem, A::AbstractOperator; samplingstyle=SamplingStyle(ap), options...) =
    solver(solverstyle, samplingstyle, ap, A; options...)
function solver(solverstyle::ProductSolverStyle, samplingstyle::ProductSamplingStyle, ap::ApproximationProblem, A::AbstractOperator; options...)
	S = sampling_operator(ap)
    solvere, sse, ape, Ae, Se = components(solverstyle), components(samplingstyle), factors(ap), components(A), factors(S)
    @assert length(solvere) == length(sse) == length(ape) == length(Ae) == length(Se)
    TensorProductOperator(
		map( (solversi, ssi, api, Ai, Si)->solver(solversi, api, Ai; S=Si, samplingstyle=ssi, options...), solvere, sse, ape, Ae, Se)...
	)
end



solver(style::AZStyle, ap::ApproximationProblem, A::AbstractOperator; options...) =
    solver(style, ap, A, AZ_Zt(ap; options...); options...)

function solver(::AZStyle, ap::ApproximationProblem, A::AbstractOperator, Zt::AbstractOperator;
            B=nothing, smallcoefficients=false, smallcoefficients_atol=NaN, smallcoefficients_rtol=NaN, verbose=false, options...)
    if smallcoefficients
        w = BasisFunctions.quadweights(sampling_grid(ap; options...), measure(ap))
        normF = abs(sqrt(sum(w .* B.^2)))
        if !isnan(smallcoefficients_rtol)
            verbose && println("Change smallcoefficients relative tolerance to absolute tolerance rtol*||f||")
            smallcoefficients_atol = smallcoefficients_rtol*normF
            smallcoefficients_rtol = NaN
        end

        AZSolver(A, Zt; smallcoefficients=smallcoefficients, smallcoefficients_rtol=smallcoefficients_rtol,
                smallcoefficients_atol=smallcoefficients_atol, verbose, options...)
    else
        AZSolver(A, Zt; B, verbose, options...)
    end
end


solver(style::AZSmoothStyle, ap::ApproximationProblem, A::AbstractOperator; options...) = solver(style, ap, A, AZ_Zt(ap; options...); options...)
solver(::AZSmoothStyle, ap::ApproximationProblem, A::AbstractOperator, Zt; options...) =
    solver(AZStyle(), ap, A, Zt; options..., weightedAZ=true, AZ_Cweight=defaultAZweight(A; options...))
