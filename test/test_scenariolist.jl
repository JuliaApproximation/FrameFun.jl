using BasisFunctions, LinearAlgebra, DomainSets, GridArrays, Test, StaticArrays, FrameFun

@testset begin
    B = (Fourier(11) → -1..1)^2
    Dom = disk(0.8)
    @test support(dictionary(∂x(random_expansion(extensionframe(B, Dom)))))≈Dom

    # @test SamplingStyle(ExtensionFramePlatform(FrameFun.ProductPlatform(FourierPlatform(),FourierPlatform()),(0.0..0.5)^2)) ==
    #     ProductSamplingStyle(OversamplingStyle(),OversamplingStyle())

    @test dictionary(_ap(Fourier(100),0.0..0.5)) == extensionframe(Fourier(100),0.0..0.5)
    # Test defaults
    @test AZ_A(Fourier(100),0.0..0.5; samplingstyle=OversamplingStyle(),oversamplingfactor=2,normalizedsampling=false)≈
        AZ_A(Fourier(100),0.0..0.5)

    ap0 = _ap(Fourier(100), 0.0..0.5; samplingstyle=OversamplingStyle(), oversamplingfactor=2, normalizedsampling=false)
    A = AZ_A(ap0)
    L = samplingparameter(ap0)
    a1 = evaluation(Fourier(L),interpolation_grid(Fourier(L)))
    a2 = IndexRestriction(GridBasis(interpolation_grid(Fourier(L))),1:200)
    a3 = BasisFunctions.FourierIndexExtension(Fourier(100),Fourier(L))
    @test A≈a2*a1*a3

    ap1 = approximationproblem(Fourier(100), 0.0..0.5; samplingstyle=OversamplingStyle(), oversamplingfactor=2, normalizedsampling=true)
    A = AZ_A(ap1)
    L = samplingparameter(ap1)
    a1 = evaluation(Fourier(L),interpolation_grid(Fourier(L)))
    a2 = IndexRestriction(GridBasis(interpolation_grid(Fourier(L))),1:200)
    a3 = BasisFunctions.FourierIndexExtension(Fourier(100),Fourier(L))
    a4 = ScalingOperator(dest(a2),1/sqrt(L))
    @test A≈a4*a2*a1*a3

    # Test differnence sampling_normalization
    ap2 = approximationproblem(Fourier(100),0.0..0.5;
            samplingstyle=OversamplingStyle(), oversamplingfactor=2, normalizedsampling=false)
    s = svdvals(AZ_A(ap2))
    N = sampling_weights(FrameFun.NormalizedSampling(OversamplingStyle()),approximationproblem(Fourier(100),0.0..0.5))
    L = samplingparameter(ap2)
    @test L == 399
    @test N[1] ≈ 1/sqrt(L)
    @test 47==sum(s.<.1*sqrt(L))
    @test 49==sum(s.>.9*sqrt(L))

    s = svdvals(AZ_A(Fourier(100),0.0..0.5;oversamplingfactor=2,normalizedsampling=true))
    @test 47==sum(s.<.1)
    @test 49==sum(1+1e-10 .> s .>.9)

    # right side should be appropriately normalized
    @test sample_data(exp, Fourier(100), 0.0..0.5; oversamplingfactor=2,normalizedsampling=true)≈
        1/sqrt(399)*GridSampling(subgrid(interpolation_grid(Fourier(399)),0.0..0.5))*exp

    # left side should be appropriately normalized
    @test AZ_A(Fourier(100),0.0..0.5;samplingstyle=OversamplingStyle(),oversamplingfactor=2,normalizedsampling=true)≈
        discretization(exp, Fourier(100), 0.0..0.5; oversamplingfactor=2,normalizedsampling=true)

    @test sample_data(exp, Fourier(100),0.0..0.5; oversamplingfactor=2,normalizedsampling=false)≈
        GridSampling(subgrid(interpolation_grid(Fourier(399)),0.0..0.5))*exp

    @test AZ_A(Fourier(100), 0.0..0.5; samplingstyle=OversamplingStyle(),oversamplingfactor=2,normalizedsampling=false)≈
        discretization(exp, Fourier(100), 0.0..0.5; oversamplingfactor=2, normalizedsampling=false)


    F = approximate(exp, Fourier(100), 0.0..0.5; samplingstyle=OversamplingStyle(),oversamplingfactor=2,normalizedsampling=true)
    @test abs(F[1](.2)-exp(.2) )< 1e-12

    F = approximate(exp, Fourier(100), 0.0..0.5; samplingstyle=OversamplingStyle(),oversamplingfactor=2,normalizedsampling=false)
    @test abs(F[1](.2)-exp(.2) )< 1e-12

    P = WeightedSumPlatform(FourierPlatform(),x->1,x->sqrt(x))
    ap = approximationproblem(P,(10,10))
    azdual(ap)
end
