module test_suite_applications
using BasisFunctions, DomainSets, FrameFun
if VERSION < v"0.7-"
    using Base.Test
else
    using Test
end


@testset "WeightedSumPlatform" begin
    # Simple platform construction test
    D = 0.0..0.5
    i = 7
    P2a = BasisFunctions.fourier_platform(oversampling=4)
    P2 = FrameFun.extension_frame_platform(P2a,D)
    WSP = FrameFun.WeightedSumPlatform([x->sqrt(x),x->1],P2)
    f = x->sqrt(x)*(1-x)-exp(x)
    AZS = AZSolver(A(WSP,i),Zt(WSP,i))
    FSP = DictFun(primal(WSP,i),AZS*f)
    rgrid = randomgrid(D,200)
    Fval = FSP(rgrid)
    # TODO: type based on space of F
    fval = sample(rgrid,f,eltype(FSP))
    abserror= sum(abs.(Fval-fval))/200
    @test abserror<1e-7
end
end
