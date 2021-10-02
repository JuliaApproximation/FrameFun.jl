
import BasisFunctions: DualType

"Trait type to select a dual dictionary suitable for the AZ algorithm."
struct AZDual <: DualType end

dual(::AZDual, dict::Dictionary, measure; options...) =
    azdual(dict, measure; options...)

"""
The dual that is used to create an AZ `Z` matrix.
"""
azdual(ss::DiscreteStyle, ap::ApproximationProblem; options...) =
    dualdictionary(platform(ap), platformparameter(ap), discretemeasure(ap); samplingstyle=ss, options...)
azdual(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    dualdictionary(platform(ap), platformparameter(ap), measure(ap); samplingstyle=ss, options...)

export azdual
@ap_interface azdual
@ap_sampling_property azdual

## The AZ algorithm
export AZ_A
@ap_interface AZ_A
AZ_A(ap::ApproximationProblem; options...) = discretization(ap; options...)


export AZ_Z
@ap_interface AZ_Z
AZ_Z(ap::ApproximationProblem; options...) = dualdiscretization(ap; options...)


export AZ_Zt
@ap_interface AZ_Zt
AZ_Zt(ap::ApproximationProblem; options...) = AZ_Z(ap; options...)'


export plungeoperator
@ap_interface plungeoperator
function plungeoperator(ap::ApproximationProblem; options...)
    A = AZ_A(ap; options...)
    Zt = AZ_Zt(ap; options...)
    I - A*Zt
end


export plungematrix, firstAZstepoperator
@ap_interface plungematrix
function plungematrix(ap::ApproximationProblem; options...)
    A = AZ_A(ap; options...)
    P = plungeoperator(ap; options...)
    P * A
end
const firstAZstepoperator = plungematrix


export plungerank
@ap_interface plungerank
function plungerank(ap::ApproximationProblem;
            rankestimate = 40,
            options...)
    C = plungematrix(ap; options...)
    Q = rSVD_solver(C; rankestimate = rankestimate, options...)
    length(Q.Sinv)
end
