
## TODO: clean up this file

include("../solvers/tridiagonalsolver.jl")
solver(::TridiagonalProlateStyle, ap, A; scaling = Zt_scaling_factor(dictionary(ap), A)
    , options...) =
    FE_TridiagonalSolver(A, scaling; options...)
# They have to do with a normalization of the
# sampling operators.
Zt_scaling_factor(S::Dictionary, A) = length(supergrid(grid(dest(A))))
Zt_scaling_factor(S::DerivedDict, A) = Zt_scaling_factor(superdict(S), A)
Zt_scaling_factor(S::ChebyshevT, A) = length(supergrid(grid(dest(A))))/2



oversampling_grid(dict::ExtensionFrame, L) = subgrid(oversampling_grid(superdict(dict),L), support(dict))

discretemeasure(ss::SamplingStyle, platform::ExtensionFramePlatform, param, ap_old; options...) =
    restrict(discretemeasure(ss, platform.basisplatform, param,
        (ap=approximationproblem(platform.basisplatform, param);setsamplingparam!(ap,samplingparameter(ap_old));ap)
            ; options...), platform.domain)

azdual_dict(sstyle::SamplingStyle, platform::ExtensionFramePlatform, param, L, measure::Measure; options...) =
   extensionframe(azdual_dict(sstyle, platform.basisplatform, param, L, supermeasure(measure); options...), platform.domain)

correct_sampling_parameter(platform::ExtensionFramePlatform, param, L_trial; options...) =
   correct_sampling_parameter(platform.basisplatform, param, L_trial; options...)

azdual_dict(sstyle::SamplingStyle, platform::AugmentationPlatform, param, L, measure::Measure; options...) =
   MultiDict([azdual_dict(sstyle, platform.basis, param, L, measure; options...), platform.functions])

oversampling_grid(samplingstyle::SamplingStyle, platform::AugmentationPlatform, param, L; options...) =
   oversampling_grid(samplingstyle, platform.basis, param, L; options...)

correct_sampling_parameter(platform::AugmentationPlatform, param, L; options...) =
    correct_sampling_parameter(platform.basis, param, L; options...)

function azdual_dict(sstyle::SamplingStyle, platform::WeightedSumPlatform, param, L, measure::Measure; options...)
    denom = (x...)->sum(map(w->abs(w(x...))^2, platform.weights))
    # TODO: discuss, what is the relation between param, L of a platform and platform.P
    MultiDict([((x...)->(platform.weights[j](x...)/denom(x...))) * azdual_dict(sstyle, platform.P, param, L, measure; options...) for j=1:length(platform.weights)])
end

oversampling_grid(samplingstyle::SamplingStyle, platform::WeightedSumPlatform, param::NTuple, L; options...) =
  oversampling_grid(samplingstyle, platform.P, param[1], L; options...)

azdual_dict(sstyle::SamplingStyle, dict::ExtensionFrame, L, measure::Measure; options...) =
   extensionframe(support(dict), azdual_dict(sstyle, superdict(dict), L, supermeasure(measure); options...),)

function azdual_dict(dict::MultiDict, measure::Measure; options...)
    dictionary = dict.dicts[1].superdict
    weights = map(weightfunction, components(dict))
    denom = (x...)->sum(map(w->abs(w(x...))^2, weights))
    MultiDict([((x...)->(weights[j](x...)/denom(x...))) * azdual_dict(dictionary, measure; options...) for j=1:length(weights)])
end

function deduce_samplingparameter(ss::OversamplingStyle, platform::WeightedSumPlatform, param::NTuple; options...)
    mix_samplingparameters(platform, map(parami->FrameFun.deduce_oversampling_parameter(ss, platform.P, parami; options...), param))
end

mix_samplingparameters(platform::WeightedSumPlatform, Ls::NTuple) =
    round.(Int,(sum(prod.(Ls))/prod(Ls[1]))^(1/length(Ls[1]))     .*Ls[1])


correct_sampling_parameter(platform::WeightedSumPlatform, param, L; options...) =
    correct_sampling_parameter(platform.P, param, L; options...)
