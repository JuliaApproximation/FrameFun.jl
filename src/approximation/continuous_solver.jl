
function continuous_approximation_operator(dest::ExtensionFrame; sampling_factor=1, solver=DirectSolver, options...)
    # since the other one is not very efficient (this one isnt either), concider this not as a general case
    (sampling_factor ≈ 1) &&
        (return ContinuousSolverPlan(solver(MixedGram(dest; options...); options...), continuous_normalization(dest; options...)))
    src = resize(dest, round(Int, sampling_factor*length(dest)))
    println(sampling_factor)
    ContinuousSolverPlan(solver(MixedGram(dest, src; options...); options...), continuous_normalization(src; options...))
end

continuous_normalization(set::Dictionary; options...) = DualGram(set; options...)
continuous_normalization(frame::ExtensionFrame; options...) = DualGram(basis(frame); options...)

struct ContinuousSolverPlan{T} <: DictionaryOperator{T}
    src                     :: Dictionary
    dest                    :: Dictionary
    mixedgramsolver         :: FE_Solver
    normalizationofb        :: DictionaryOperator

    scratch                 :: Vector{T}
    mixedgram               :: DictionaryOperator
    ContinuousSolverPlan{T}(src::Dictionary, dest::Dictionary, mixedgramsolver::FE_Solver, normalizationofb::DictionaryOperator) where {T} =
        new(src, dest, mixedgramsolver, normalizationofb, zeros(T, length(src)), operator(mixedgramsolver))
end

ContinuousSolverPlan(solver::FE_Solver{T}, normalization::DictionaryOperator) where {T} =
    ContinuousSolverPlan(src(solver), dest(solver), solver, normalization)

ContinuousSolverPlan(src::Dictionary, dest::Dictionary, solver::FE_Solver{T}, normalization::DictionaryOperator) where {T} =
    ContinuousSolverPlan{T}(src, dest, solver, normalization)

function apply!(s::ContinuousSolverPlan, coef_dest, coef_src)
    apply!(s.normalizationofb, s.scratch, coef_src)
    coef_dest[:] = s.mixedgramsolver*s.scratch
end
