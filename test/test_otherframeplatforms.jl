
using BasisFunctions, FrameFun.AugmentationPlatforms, FrameFun.WeightedSumPlatforms, Test, FrameFun.Platforms


@testset  "AugmentationPlatform" begin
    plat1 = AugmentationPlatform(platform(Fourier(10)),Fourier(5))
    plat2 = ONB_plus_K(platform(Fourier(10)),sqrt, 5)
    plat3 = AugmentationPlatform(platform(ChebyshevT(10)),ChebyshevT(5))
    plat4 = ONB_plus_K(platform(ChebyshevT(10)),sqrt, 5)
    plat5 = AugmentationPlatform(platform(ChebyshevT(10)),Fourier(5))
    plat6 = AugmentationPlatform(platform(Fourier(10)),ChebyshevT(5))

    for plat in (plat1,plat2,plat3,plat4,plat5,plat6)
        @test SamplingStyle(plat) == OversamplingStyle()
        @test DictionaryStyle(plat) == FrameStyle()
        @test SolverStyle(plat,OversamplingStyle()) == AZStyle()
    end

    for plat in (plat1,plat2,plat6)
        @test measure(plat) == FourierMeasure()
    end
    for plat in (plat3,plat4,plat5)
        @test measure(plat) == ChebyshevMeasure()
    end

    dict1 = dictionary(plat1,20)
    dict2 = dictionary(plat2,20)
    dict3 = dictionary(plat3,20)
    dict4 = dictionary(plat4,20)
    dict5 = dictionary(plat5,20)
    dict6 = dictionary(plat6,20)

    for dict in (dict1,dict2,dict3,dict4,dict5,dict6)
        @test dict isa MultiDict
        @test dimensions(dict) == [20,5]
    end

    @test_throws Exception dualdictionary(plat1,20,measure(plat1))
    @test_throws Exception dualdictionary(plat2,20,measure(plat2))
    @test_throws Exception dualdictionary(plat3,20,measure(plat3))
    @test_throws Exception dualdictionary(plat4,20,measure(plat4))
    @test_throws Exception dualdictionary(plat5,20,measure(plat5))
    @test_throws Exception dualdictionary(plat6,20,measure(plat6))

    # TODO Product platform tests

end

@testset "WeightedSumPlatforms" begin
    plat = WeightedSumPlatform(platform(Fourier(10)),x->1., sqrt)
    @test SamplingStyle(plat) == OversamplingStyle()
    @test DictionaryStyle(plat) == FrameStyle()
    @test SolverStyle(plat,OversamplingStyle()) == AZStyle()
    @test measure(plat) == FourierMeasure()
    dict = dictionary(plat,(10,10))
    @test dict isa MultiDict
    @test dimensions(dict) == [10,10]
    ddict =  dualdictionary(plat,(10,10),FourierMeasure())
    @test ddict isa MultiDict
    @test dimensions(ddict) == [10,10]
    @test param_first(plat) == (10,10)

end
