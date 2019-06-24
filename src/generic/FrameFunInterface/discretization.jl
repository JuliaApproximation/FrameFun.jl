discretization(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    discretization(ss, ap, haskey(options,:S) ? options[:S] : samplingoperator(ss, ap; options...); options...)
discretization(::SamplingStyle, ap::ApproximationProblem, S::AbstractOperator; options...) =
    apply(S, dictionary(ap); options...)
discretization(ss::ProductSamplingStyle, ap::ApproximationProblem, S::AbstractOperator; options...) =
    tensorproduct( elements(discretization, ss, ap,
        BasisFunctions.productelements(S); options...)... )

# The discretization routine can also take f as an argument
discretization(f, dict::Dictionary, args...; options...) =
    discretization(f, approximationproblem(dict, args...); options...)
discretization(f, platform::Platform, args...; options...) =
    discretization(f, approximationproblem(platform, args...); options...)
discretization(f, ap::ApproximationProblem; samplingstyle=SamplingStyle(ap), options...) =
    discretization(f, samplingstyle, ap; options...)
discretization(f, ss::SamplingStyle, ap::ApproximationProblem; options...) =
    discretization(f, ss, ap, haskey(options,:S) ? options[:S] : samplingoperator(ss, ap; options...); options...)
function discretization(f, ss::SamplingStyle, ap::ApproximationProblem, S::AbstractOperator; options...)
    A = discretization(ss, ap, S; options...)
    B = apply(S, f; options...)
    A, B
end

function dualdiscretization(ss::SamplingStyle, ap::ApproximationProblem; options...)
    dualdict = haskey(options, :dualdict) ? options[:dualdict] : azdual_dict(ss, ap; options...)
    S = haskey(options, :S) ? options[:S] : dualsamplingoperator(ss, ap; dualdict=dualdict, options...)
    dualdiscretization(ss, ap, S, dualdict; options...)
end

function dualdiscretization(ss::SamplingStyle, ap::ApproximationProblem, S::AbstractOperator; options...)
    dualdict = haskey(options, :dualdict) ? options[:dualdict] : azdual_dict(ss, ap; options...)
    dualdiscretization(ss, ap, S, dualdict; options...)
end

dualdiscretization(ss::SamplingStyle, ap::ApproximationProblem, Stilde::AbstractOperator, dualdict::Dictionary; options...) =
    apply(Stilde, dualdict; options...)

dualdiscretization(ss::ProductSamplingStyle, ap::ApproximationProblem, S::AbstractOperator, dualdict::TensorProductDict; options...) =
    tensorproduct( elements(dualdiscretization, ss, ap,
        BasisFunctions.productelements(S), elements(dualdict); options...)... )
