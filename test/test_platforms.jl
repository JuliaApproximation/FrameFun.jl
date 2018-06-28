module test_suite_applications
using BasisFunctions
using FrameFun
using FrameFun: FrameFun.spline_util_restriction_operators, boundary_support_grid, relative_indices, restriction_operator
using Base.Test
using StaticArrays
using Domains


@testset "WeightedSumPlatform" begin
    # Simple platform construction test
    D = interval(0,0.5)
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
    @test abserror<1e-8
end
end
    
