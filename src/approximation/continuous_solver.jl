
function continuous_approximation_operator(set::FunctionSet; sampling_factor=1, solver=FrameFun.FE_DirectSolver, options...)
    # since the other one is not very efficient (this one isnt either), concider this not as a general case
    (sampling_factor ≈ 1) &&
        (return ContinuousSolverPlan(set, solver(MixedGram(set; options...)), continuous_normalization(set; options...)))
    set2 = resize(set, round(Int, sampling_factor*length(set)))

    ContinuousSolverPlan(set, solver(MixedGram(set; options...)), continuous_normalization(set; options...)
end

continuous_normalization(set::FunctionSet; options...) = DualGram(set, options...)
continuous_normalization(frame::ExtensionFrame; options...) = DualGram(set, options...)

immutable ContinuousSolverPlan{T} <: AbstractOperator{T}
    src                     :: FunctionSet
    dest                    :: FunctionSet
    mixedgramsolver         :: FE_Solver
    normalizationofb        :: AbstractOperator
    scratch                 :: Array{T,1}
    ContinuousDirectSolver{ELT}(src::FunctionSet, dest::FunctionSet, mixedgramsolver::FE_Solver, normalizationofb::AbstractOperator) where {ELT,SOLVER<:FE_Solver} =
        new(set, mixedgramsolver, normalizationofb, zeros(set))
end

function apply!(s::ContinuousSolverPlan, coef_dest, coef_src)
    apply!(s.normalizationofb, s.scratch, coef_src)
    coef_dest[:] = s.mixedgramsolver*s.scratch
end
